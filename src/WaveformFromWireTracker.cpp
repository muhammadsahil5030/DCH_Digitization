#include "WaveformFromWireTracker.h"

#include "DCHUtilities.h"
#include "DCHFFTNoise.h"
#include "DCHXT2DLUT.h"

#include "TMath.h"
#include "TF1.h"

// STL
#include <algorithm>
#include <sstream>
#include <cmath>
#include <unordered_map>

// Gaudi
#include <GaudiKernel/MsgStream.h>

#include "DD4hep/Detector.h"
#include "DD4hep/DetElement.h"

///////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// 	 constructor       ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
// -- KeyValues("name of the variable that holds the name of the collection exposed in the python steering file",
// {"default name for the collection"}),
WaveformFromWireTracker::WaveformFromWireTracker(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc,
		{
		    //input collection
            KeyValues("InputClusterCollection", {"GasIonizationClusters"}),
            KeyValues("HeaderName", {"EventHeader"}),
        },
		{
			// hWaveformDigiL
			KeyValues("DigiWaveformLCollection",   {"SenseWire_DigiWaveformsL"}),
			// hWaveformDigiR
			KeyValues("DigiWaveformRCollection",   {"SenseWire_DigiWaveformsR"})
		}
	) 
	{
		m_geoSvc = serviceLocator()->service(m_geoSvcName);	// Access detector geometry service
		m_uidSvc = serviceLocator()->service(m_uidSvcName);	// Access reproducible unique-ID seed service
	}

// Random number generator
TRandom3 WaveformFromWireTracker::CreateRandomEngine(
		const edm4hep::EventHeaderCollection& headers) const {
  const auto seed = m_uidSvc->getUniqueID(headers, this->name());
  TRandom3 rng(seed);
  return rng;
}

//==================================================================================//
// 						Function Sampling avalache charge							//
//==================================================================================//
double WaveformFromWireTracker::avalancheCharge(TRandom3& myRandom) const {
	if (!m_polya) {
		return 0.0;
	}
	const double q = m_polya->GetRandom(&myRandom);
	return (q > 0.0) ? q : 0.0;
}

//=================================================================================//
//  				Fuction of Single-electron pulse shape						   //
//=================================================================================//
// t_ns   : observation time (ns)
// t0_ns  : electron arrival time (ns)
// q      : avalanche charge (arbitrary units, e.g. from Polya)
// Shape: double exponential CR-RC-like pulse
double WaveformFromWireTracker::singleElectronPulse(double t_ns, double t0_ns, double q) const {
	const double dt = t_ns - t0_ns;
	// No signal before the electron arrival
	if (dt <= 0.0) {
		return 0.0;
	}
 
	const double tau_r = m_pulseRiseTime_ns.value();
	const double tau_f = m_pulseFallTime_ns.value();
	// Safety: if parameters are not set, just return a delta-like pulse
	if (tau_r <= 0.0 || tau_f <= 0.0 || tau_r == tau_f) {
		//return q;
		return 0.0;
	}

	// Double exponential shape
	const double term_f = std::exp(-dt / tau_f);
	const double term_r = std::exp(-dt / tau_r);
	// Optional normalization factor so peak ~ q (roughly)
	const double norm = 1.0 / (tau_f - tau_r);
	const double q_eff = q * m_pulseAmplitudeScale.value();

	return q_eff * norm * (term_f - term_r);
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////       initialize       ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
StatusCode WaveformFromWireTracker::initialize() 
{
	// Check for UniqueIDGen and Geometery service
	if (!m_uidSvc) {
		error() << "Unable to get UniqueIDGenSvc" << endmsg;
		return StatusCode::FAILURE;
	}
	if (!m_geoSvc) {
		error() << "Unable to get GeoSvc" << endmsg;
		return StatusCode::FAILURE;
	}

	// Retrieve the subdetector
	std::string detectorName(m_detectorName.value());
	if (0 == m_geoSvc->getDetector()->detectors().count(detectorName)) {
		error() << "Detector <<" << detectorName <<  ">> does not exist." << endmsg;
		return StatusCode::FAILURE;
	}

	dd4hep::DetElement detectorElement = m_geoSvc->getDetector()->detectors().at(detectorName);

	// Retrieve the wire-tracker extension using the DD4hep base extension key.
	auto* wireTracker_info = detectorElement.extension<dd4hep::rec::WireTracker_info_struct>();
	if (!wireTracker_info) {
		error() << "No WireTracker_info_struct extension was attached to detector <<"
			<< detectorName << ">>." << endmsg;
		return StatusCode::FAILURE;
	}

	// Cast the generic wire-tracker extension to the concrete DCH geometry helper.
	dch_data = dynamic_cast<dd4hep::rec::DCH_info*>(wireTracker_info);
	if (!dch_data || !dch_data->IsValid()) {
		error() << "No valid wire-tracker geometry extension was found for detector <<"
				<< detectorName << ">>." << endmsg;
		return StatusCode::FAILURE;
	}

	// Retrieve the readout associated with the detector element (subdetector)
	dd4hep::SensitiveDetector sensitiveDetector = m_geoSvc->getDetector()->sensitiveDetector(detectorName);
	if (!sensitiveDetector.isValid()) {
		error() << "No valid Sensitive Detector was found for detector <<" << detectorName << ">>." 
		<< endmsg;
		return StatusCode::FAILURE;
	}

	dd4hep::Readout readout = sensitiveDetector.readout();
	// set the cellID decoder
	m_decoder = readout.idSpec().decoder();
	if (!m_decoder) {
		error() << "Failed to retrieve CellID decoder from detector readout <<" << readout.name() << ">>." 
		<< endmsg;
		return StatusCode::FAILURE;
	}

	// Read DCH half-length from XML constant (DCH_standalone_o1_v02.xml)
	dd4hep::Detector& det = dd4hep::Detector::getInstance();

	if (!det.constantAsDouble("DCH_gas_Lhalf")) {
		error() << "Missing DD4hep constant: DCH_gas_Lhalf" << endmsg;
		return StatusCode::FAILURE;
	}

	m_halfChamberLength_mm = det.constantAsDouble("DCH_gas_Lhalf") / MM_TO_CM; // returns in mm
	info() << "[WaveformFromWireTracker] m_halfChamberLength_mm = " << m_halfChamberLength_mm << " mm" 
	<< endmsg;

	//=================================================================================//
	//							  initialize x-t relation  							 //
	//=================================================================================// 
	// Load the x-t lookup table
	try {
		// Create xt LUT object
		m_xt2d = std::make_unique<DCHXT2DLUT>(m_xtFileName.value(), "xt_mean", "xt_error");
	} catch (const std::exception& e) {
		error() << e.what() << endmsg;
		return StatusCode::FAILURE;
	}
	info() << "[XT2D] Loaded XT LUT from " << m_xtFileName.value() << endmsg;

	//=================================================================================//
	//						    Load FFT-noise template    							 //
	//=================================================================================//
  	if (m_enableFFTNoise.value()) {
		std::string err;
		// Load FFT magnitude template + normalization + frequency
		const bool noisefile = dchfft::loadFFTNoiseTemplate(
			m_noiseFileName.value(),	// ROOT file containing the FFT-noise template
			m_noiseDirName.value(),	
			m_fftMag, 	//Magnitude
			m_fftAmpNorm,
			m_fftMaxFreq,
			m_fftSize,
			&err
		);

		if (!noisefile) {
			error() << "[FFTNoise] " << err << endmsg;
			return StatusCode::FAILURE;
		}
		// Mark noise system ready for waveform generation
		m_fftNoiseReady = true;
		info() << "[FFTNoise] Loaded fft_noise:"
		<< " size=" << m_fftSize
		<< " maxFreq=" << m_fftMaxFreq
		<< " ampNorm=" << m_fftAmpNorm
		<< endmsg;  
	}

	//=================================================================================//
	//							  initialize Polya gain  							   //
	//=================================================================================//
	// Polya PDF for gas gain:
	// P(Q) ~ ((1+theta)/(G * Gamma(theta+1))) *
	//        [Q*(1+theta)/G]^theta * exp( - Q*(1+theta)/G )
	// Parameters:
	//   [0] normalization (kept at 1; ROOT will renormalize as needed)
	//   [1] mean gas gain G
	//   [2] theta (Polya shape parameter)
	// Range (1e3, 1e7) is taken from the original GMCT implementation
	// and can be tuned later if needed.
	//
	m_polya = new TF1("Polya",
		"((1.+[2])/([1]*TMath::Gamma([2]+1))) * pow(x*(1.+[2])/[1],[2]) * exp(-x*(1.+[2])/[1])",
		1.e3, 1.e7);
	m_polya->SetParameter(0, 1.0);
	m_polya->SetParameter(1, m_GasGain.value());
	m_polya->SetParameter(2, m_PolyaTheta.value());
	//-----------------------------------------------------------------------------------//

	std::stringstream ss;
	PrintConfiguration(ss);
	info() << ss.str().c_str() << endmsg;

	// Create histograms for debugging:
	if (m_create_debug_histos.value()) {
		m_histSvc = service("THistSvc");
		if (!m_histSvc) {
			error() << "Unable to locate Histogram Service" << endmsg;
			return StatusCode::FAILURE;
		}

		hAvalancheQ = new TH1F("hAvalancheQ", "Avalanche per Charge; Q  (a.u.); Entries", 100, 0, 150.e4);
		if (m_histSvc->regHist("/rec/hAvalancheQ", hAvalancheQ).isFailure()) {
			error() << "Couldn't register hAvalancheQ histogram" << endmsg;
		}

		const double pulseTmin = m_PulseStart_ns.value();
		const double pulseTmax = m_PulseWindow_ns.value();
		const double dtPulse_ns = m_pulseBinSize_ns.value();
		const int nBinsSignalPulse = std::max(
			1, static_cast<int>(std::lround((pulseTmax - pulseTmin) / dtPulse_ns)));

		hSignalPulse = new TH1F("hSignalPulse",
								"Single-electron pulse; t [ns]; amplitude [V]",
								nBinsSignalPulse, pulseTmin, pulseTmax);
		if (m_histSvc->regHist("/rec/hSignalPulse", hSignalPulse).isFailure()) {
			error() << "Couldn't register hSignalPulse histogram" << endmsg;
		}

		const double wfTmin = m_waveformStartTime_ns.value();
		const double wfTmax = wfTmin + m_waveformTimeWindow_ns.value();
		
		const double dtAnalog_ns = m_waveformBinSize_ns.value();
		const double dtDigi_ns = m_adcSamplePeriod_ns.value();


		const int nBinsAnalogWf  = static_cast<int>(std::lround((wfTmax - wfTmin) / dtAnalog_ns));
		const int nBinsDigiWf = static_cast<int>(std::lround((wfTmax - wfTmin)/dtDigi_ns));

		hWaveform = new TH1F("hWaveform", "Analog Waveform; t [ns]; Amplitude", nBinsAnalogWf, wfTmin, wfTmax);
		if (m_histSvc->regHist("/rec/hWaveform", hWaveform).isFailure()) {
			error() << "Couldn't register hWaveform histogram" << endmsg;
		}
		hWaveformL = new TH1F("hWaveformL", "Waveform at wire end 1; t [ns]; Amplitude", nBinsAnalogWf, wfTmin, wfTmax);
		if (m_histSvc->regHist("/rec/hWaveformL", hWaveformL).isFailure()) {
			error() << "Couldn't register hWaveformL histogram" << endmsg;
		}
		hWaveformR = new TH1F("hWaveformR","Waveform at wire end 2; t [ns]; Amplitude",nBinsAnalogWf, wfTmin, wfTmax);
		if (m_histSvc->regHist("/rec/hWaveformR", hWaveformR).isFailure()) {
			error() << "Couldn't register hWaveformR histogram" << endmsg;
		}
		//-------------------------------------------------------------------------------------------//

		hWaveformElecL = new TH1F("hWaveformElecL", "Waveform after electronics (end 1); t [ns]; Amplitude",
				nBinsAnalogWf, wfTmin, wfTmax);
		if (m_histSvc->regHist("/rec/hWaveformElecL", hWaveformElecL).isFailure()) {
			error() << "Couldn't register hWaveformElecL histogram" << endmsg;
		}
		hWaveformElecR = new TH1F("hWaveformElecR","Waveform after electronics (end 2); t [ns]; Amplitude",
								nBinsAnalogWf, wfTmin, wfTmax);
		if (m_histSvc->regHist("/rec/hWaveformElecR", hWaveformElecR).isFailure()) {
			error() << "Couldn't register hWaveformElecR histogram" << endmsg;
		}
		//------------------------------------------------------------------------------------------//

		// Debug histogram for FFT noise (single example waveform)
		hNoise = new TH1F("hNoise", "FFT Noise; t [ns]; Amplitude",
				nBinsAnalogWf, wfTmin, wfTmax);
		if (m_histSvc->regHist("/rec/hNoise", hNoise).isFailure()) {
			error() << "Couldn't register hNoise histogram" << endmsg;
		}

		// Digitized (reconstructed) waveform in mV
		hWaveformDigiL = new TH1F("hWaveformDigiL","Digitized waveform at end 1; t [ns]; voltage [mV]",
								nBinsDigiWf, wfTmin, wfTmax);
		if (m_histSvc->regHist("/rec/hWaveformDigiL", hWaveformDigiL).isFailure()) {
			error() << "Couldn't register hWaveformDigiL histogram" << endmsg;
		}
		hWaveformDigiR = new TH1F("hWaveformDigiR","Digitized waveform at end 2; t [ns]; voltage [mV]",
								nBinsDigiWf, wfTmin, wfTmax);
		if (m_histSvc->regHist("/rec/hWaveformDigiR", hWaveformDigiR).isFailure()) {
			error() << "Couldn't register hWaveformDigiR histogram" << endmsg;
		}
	} //End of histogram definition

	return StatusCode::SUCCESS;
} //end of Initialize

///////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////       	operator()       //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
std::tuple< edm4hep::RawTimeSeriesCollection, edm4hep::RawTimeSeriesCollection>
WaveformFromWireTracker::operator()(
	const edm4hep::SimTrackerHitCollection& input_clusters,
	const edm4hep::EventHeaderCollection&   headers
) const 
{

	// initialize random engine
	auto myRandom = this->CreateRandomEngine(headers);

	// Create the collections we are going to return
	edm4hep::RawTimeSeriesCollection out_Digitized_waveformsL;
	edm4hep::RawTimeSeriesCollection out_Digitized_waveformsR;

	//Varible initailizations:
	int ilayer = -1;
	int nphi = -1;

	// Local analog waveform definition used for ALL cells
	const double wfTmin = m_waveformStartTime_ns.value();
	const double wfTmax = wfTmin + m_waveformTimeWindow_ns.value();
	const double waveformBinSize_ns = m_waveformBinSize_ns.value();
	
	const int nWaveformBins = std::max(
		1, static_cast<int>(std::lround((wfTmax - wfTmin) / waveformBinSize_ns)));

	std::vector<double> waveformTimeAxis(nWaveformBins, 0.0);
	for (int i = 0; i < nWaveformBins; ++i) {
		waveformTimeAxis[i] = wfTmin + (i + 0.5) * waveformBinSize_ns;
	}

  	// Container used to accumulate all electron-level information for one readout cell.
	struct CellAccum {
		std::vector<double> times_ns;   // Absolute electron arrival times [ns]
		std::vector<double> charges;    // Avalanche charge per electron
		std::vector<double> wireZ_mm;   // Electron position along the wire [mm]

		float wireStereoAngle = 0.f;	// Representative stereo angle of the sense wire for this cell
	};
	// Map from cellID to the accumulated electron information of that cell.
  	std::unordered_map<uint64_t, CellAccum> cellAccumMap;

  	// Each input object represents one generated cluster
	// loop over hit collection
	int loop_index = 0;
  	for (const auto& cluster : input_clusters) 
  	{
		loop_index++;

		// Extract basic cluster information: cell ID, position, and time.
		dd4hep::DDSegmentation::CellID cellid = cluster.getCellID();
		//auto hit_position = Convert_EDM4hepVector_to_TVector3(cluster.getPosition(), MM_TO_CM);
		const double clusterTime_ns = cluster.getTime();

		// Decode the cluster cell ID into geometry indices.
		const int superlayer = this->CalculateSuperLayerFromCellID(cellid);
		const int sector = this->CalculateSectorFromCellID(cellid);
		ilayer = this->CalculateLayerFromCellID(cellid);
		nphi = this->CalculateNphiFromCellID(cellid);

		auto hit_position = Convert_EDM4hepVector_to_XYZVector(cluster.getPosition(), MM_TO_CM);

		//Print some variable
		debug() << "Cluster #" << loop_index
			<< " layer=" << ilayer
			<< " cell=" << nphi
			<< endmsg;

		// Sample the number of electrons associated with this cluster.
		// The cluster size is sample from He experimental cluster size distribution. 
		const int cluster_size = sample_cluster_size(myRandom);
		if (cluster_size <= 0) {
			continue;
		}
    
		// Compute the local sense-wire direction. 
		// This is later used to build the transverse plane for the x-t lookup.
		//TVector3 wire_direction_ez = this->dch_data->Calculate_wire_vector_ez(ilayer, nphi);
		auto wire_direction_ez = this->dch_data->Calculate_wire_vector_ez(superlayer, ilayer, sector, nphi);
		float WireStereoAngle = 0;
		{
		const auto& layer = this->dch_data->database.at(ilayer);
		// radial middle point of the cell at Z=0
		auto cell_rave_z0 = 0.5 * (layer.radius_fdw_z0 + layer.radius_fuw_z0);
		// when building the twisted tube, the twist angle is defined as:
		//     cell_twistangle    = layer.StereoSign() * DCH_i->twist_angle
		// which forces the stereoangle of the wire to have the oposite sign
		WireStereoAngle = static_cast<float>(
			-1.0 * dch_data->StereoSign(layer) * dch_data->stereoangle_z0(cell_rave_z0)
		);
		}

		//==================================================================================//
		//                    Drift time for one already-generated cluster                  //
		//==================================================================================//
		// Store one drift time per electron belonging to this cluster.
		std::vector<double> electron_times_ns;
		electron_times_ns.reserve(cluster_size);

		// Store the projected z position of each electron along the sense wire.
		std::vector<double> electron_wireZ_mm;
		electron_wireZ_mm.reserve(cluster_size);

		auto cluster_pos = hit_position;

		auto cl_to_wire = this->dch_data->Calculate_hitpos_to_wire_vector(
			superlayer, ilayer, sector, nphi, cluster_pos);

		auto cl_on_wire = cluster_pos + cl_to_wire;

		const double z_wire_cl_mm = cl_on_wire.z() / MM_TO_CM;

		auto u = wire_direction_ez.Unit();

		ROOT::Math::XYZVector e1(1.0, 0.0, 0.0);
		e1 -= e1.Dot(u) * u;

		if (std::sqrt(e1.Dot(e1)) < 1e-9) {
		e1 = ROOT::Math::XYZVector(0.0, 1.0, 0.0);
		e1 -= e1.Dot(u) * u;
		}

		e1 = e1.Unit();

		auto e2 = u.Cross(e1).Unit();

		auto v_perp = cl_to_wire - cl_to_wire.Dot(u) * u;

		const double x_cluster_cm = v_perp.Dot(e1);
		const double y_cluster_cm = v_perp.Dot(e2);

		// Sample one drift time per electron in the cluster
		for (int ie = 0; ie < cluster_size; ++ie) {
			const double t_ns = m_xt2d->sampleTimeNs(x_cluster_cm, y_cluster_cm, myRandom);

			electron_times_ns.push_back(t_ns);
			electron_wireZ_mm.push_back(z_wire_cl_mm);
		}

		//==========================================================================================
		//                               Apply Polya distribution                                 //
		//==========================================================================================
		// Store avalanche charge for each drifted electron
		std::vector<double> electron_charges;
		electron_charges.reserve(electron_times_ns.size());
		double Qtot = 0.0;

		//for each drifted electron, sample the avalunche uisng
		//Polya distribution.
		for (size_t i = 0; i<electron_times_ns.size(); ++i)
		{
			double q = avalancheCharge(myRandom);

			electron_charges.push_back(q);
			Qtot +=q;

			if (m_create_debug_histos.value() && hAvalancheQ) {
				hAvalancheQ->Fill(q);
			}
		}

		//=======================================================================================
		// 						Accumulate electrons per cellID				       			   //
		//======================================================================================= 			
		auto& cellAccum = cellAccumMap[static_cast<uint64_t>(cellid)];

		// Store one representative stereo angle for this cell.
		cellAccum.wireStereoAngle = WireStereoAngle;

		// Store absolute electron times in the event time frame
		for (double tdrift_ns : electron_times_ns) {
			cellAccum.times_ns.push_back(clusterTime_ns + tdrift_ns);
		}

		cellAccum.charges.insert(cellAccum.charges.end(),
				electron_charges.begin(),
				electron_charges.end());

		cellAccum.wireZ_mm.insert(
			cellAccum.wireZ_mm.end(),
			electron_wireZ_mm.begin(),
			electron_wireZ_mm.end());
	
  	} // end loop over hit collection
	//**************************************************************************//
  
	// look for the first cell to save only one debug histogram 
  	uint64_t debugCellID = 0;

  	// Retrieve accumulated data per cell (all electrons)
  	for (const auto& kv : cellAccumMap) 
  	{
		const auto& acc = kv.second;
				
		if (acc.times_ns.empty()) continue;  

		// Store the electrons charge, drift time, and 
		// distance along the wire per cell:
		const auto& electron_times_ns = acc.times_ns;	   
		const auto& electron_charges  = acc.charges;	    
		const auto& electron_wireZ_mm = acc.wireZ_mm;

			
		// Local waveform buffers for THIS cell only
		std::vector<double> wfCenter(nWaveformBins, 0.0);
		std::vector<double> wfLeft(nWaveformBins, 0.0);
		std::vector<double> wfRight(nWaveformBins, 0.0);
		std::vector<double> wfElecL(nWaveformBins, 0.0);
		std::vector<double> wfElecR(nWaveformBins, 0.0);

		// Wire stereo angle
		const float WireStereoAngle = acc.wireStereoAngle;
		
		// Earliest electron arrival time in cell
		auto it_min = std::min_element(acc.times_ns.begin(), acc.times_ns.end());
		const double t_hit_ns =  *it_min;

		// Select first valid cell to store debug waveform histograms
		if (m_create_debug_histos.value() && !m_waveformFilled) {
				debugCellID = kv.first;
		}

		const bool isDebugCell = m_create_debug_histos.value() &&
			!m_waveformFilled &&
			(kv.first == debugCellID);
			

		//==========================================================================================
		//               Build analog waveform by summing the individual Pulses
		//==========================================================================================

		if (!electron_times_ns.empty()) 
	  	{
    		if (isDebugCell && hSignalPulse &&
					electron_times_ns.size() == electron_charges.size()) {

				hSignalPulse->Reset();

				// Find earliest electron arrival time in this cell
				auto it_min_e = std::min_element(electron_times_ns.begin(), 
						electron_times_ns.end());

				// Get index of the earliest electron
				const size_t idx_min = std::distance(electron_times_ns.begin(), 
						it_min_e);

				// Use charge of the earliest electron for the reference single-electron pulse
				const double q = electron_charges[idx_min];
				const double t0_ns = electron_times_ns[idx_min];

				// Fill the debug pulse histogram bin by bin in time
				for (int ibin = 1; ibin <= hSignalPulse->GetNbinsX(); ++ibin) {
						const double t = hSignalPulse->GetBinCenter(ibin);
						const double v = singleElectronPulse(t, t0_ns, q);
						hSignalPulse->SetBinContent(ibin, v);
				}
	      	}

			// Loop over all time bins of the analog waveform
			for (int ibin = 0; ibin < nWaveformBins; ++ibin) {
				const double t = waveformTimeAxis[ibin];
				double V_t = 0.0;

				// Sum contribution from every electron arriving
				for (size_t i = 0; i < electron_times_ns.size(); ++i) {
					const double t0_ns = electron_times_ns[i];
					const double q     = electron_charges[i];

					// Add this electron pulse contribution to the waveform sample
					V_t += singleElectronPulse(t, t0_ns, q);
				}
				wfCenter[ibin] = V_t;
			}

			//=========================================================================================
			//                   propagate waveform to both ends (delay + attenuation)
			//=========================================================================================
			// Stereo angle
			const double cosStereo = std::cos(WireStereoAngle);
			const double cosangle  = (std::abs(cosStereo) > 1e-6) ? std::abs(cosStereo) : 1.0;

			// Signal propagation speed along the sense wire.
			const double v_mm_per_ns = m_signalSpeed_mm_per_ns.value();
			
			// Attenuation length for signal transport along the wire.
			const double lambda_mm   = m_wireAttenuationLength_mm.value();
			const bool useAttenuation = (lambda_mm > 0.0);

			// Sum the propagated signal seen at both wire ends time-bin by time-bin.
	      	for (int ibin = 0; ibin < nWaveformBins; ++ibin) {
		  	  	const double t = waveformTimeAxis[ibin];

				double VL = 0.0;
				double VR = 0.0;

				// Each electron reaches the two ends with different delay and attenuation.
				for (size_t i = 0; i < electron_times_ns.size(); ++i) {
					// Electron production point along the wire.
						const double z_wire_mm = electron_wireZ_mm[i];

					// Propagation distance from the electron point to the left-right wire end.
						double distR_mm = (m_halfChamberLength_mm - z_wire_mm) / cosangle;
						double distL_mm = (m_halfChamberLength_mm + z_wire_mm) / cosangle;
						if (distR_mm < 0.0) distR_mm = 0.0;
						if (distL_mm < 0.0) distL_mm = 0.0;

					// Convert wire-travel distance into propagation time.
						const double tpropR_ns = (v_mm_per_ns > 0.0) ? distR_mm / v_mm_per_ns : 0.0;
						const double tpropL_ns = (v_mm_per_ns > 0.0) ? distL_mm / v_mm_per_ns : 0.0;
					
					// Left, Right-end amplitude loss from wire attenuation.
					const double attR = useAttenuation ? std::exp(-distR_mm / lambda_mm) : 1.0;
					const double attL = useAttenuation ? std::exp(-distL_mm / lambda_mm) : 1.0;

					// Local electron time inside the waveform window.
					const double t0_ns = electron_times_ns[i];

					// Left, Right-end signal = delayed and attenuated single-electron pulse.
					VR += singleElectronPulse(t, t0_ns + tpropR_ns,
							electron_charges[i] * attR);

					VL += singleElectronPulse(t, t0_ns + tpropL_ns,
							electron_charges[i] * attL);
				}

				// Store propagated waveform amplitudes at the two readout ends.
				wfLeft[ibin]  = VL;
				wfRight[ibin] = VR;
			}
	  
		}
		// ==========================================================================================//
		//                          impedance mismatch (reflection/transmission)                     //
		//===========================================================================================//
		//This function return the electronics transfer function impulse shape
		//by using the rise and fall time
		auto gmct_signalShape2 = [&](double t0_ns, double amp) -> double 
		{
			const double tauFall1 = m_elecTauFall1_ns.value();
			const double tauRise  = m_elecTauRise_ns.value();
			const double mix      = m_elecMixFraction.value();
			const double tauFall2 = m_elecTauFall2_ns.value();

			if (t0_ns <= 0.0) return 0.0;
			if (tauFall1 <= 0.0 || tauFall2 <= 0.0 || tauRise <= 0.0) return 0.0;

			// Double-exponential decay describing front-end shaping tail.
			double sign = mix * std::exp(-t0_ns / tauFall1) / tauFall1;
			sign       += (1.0 - mix) * std::exp(-t0_ns / tauFall2) / tauFall2;
			// Rising edge shaping
			sign       *= (1.0 - std::exp(-t0_ns / tauRise));

			// Normalization factor
			sign       *= (tauFall1 + tauRise) * (tauRise + tauFall2);
			sign       /= (tauFall1 * tauFall2 + tauRise * tauFall2
					+ tauFall1 * tauRise * mix - tauRise * mix * tauFall2);

			// Return shaped signal scaled by input amplitude.
			return sign * amp;

		};

		// Apply impedance mismatch + electronics shaping via convolution.
		auto apply_step10_electronics_vec = [&](const std::vector<double>& in, std::vector<double>& out)
		{
			// Initialize output waveform
			out.assign(in.size(), 0.0);
			const int nBins = static_cast<int>(in.size());
			if (nBins <= 0) return;

			const double dt_ns = waveformBinSize_ns;
					
			// Start with front-end gain scaling.
			double scale = m_frontEndGain.value();

			// Impedance mismatch modifies transmitted signal amplitude
			// (reflection losses).
			if (m_enableImpedanceMismatch.value()) {
				const double R = m_matchingRes_Ohm.value();
				const double Z = m_tubeImpedance_Ohm.value();
				const double denom = (R + Z);

				// Reflection coefficient (R-Z)/(R+Z)
				if (std::abs(denom) > 0.0) {
					const double reflect = (R - Z) / denom;
					scale *= (1.0 - reflect);
				}
			}

			// Apply global gain + transmission scaling to waveform.
			std::vector<double> x(nBins, 0.0);
			for (int i = 0; i < nBins; ++i) {
				x[i] = in[i] * scale;
			}
			// If TF disabled, just copy scaled waveform
			if (!m_enableElectronicsTF.value()) {
				out = x;
				return;
			}

			// Build finite-length impulse response kernel of electronics.
			const int nKernel = std::min( nBins, std::max(1, 
						int(std::ceil(m_elecKernelLength_ns.value() / dt_ns)) + 1));
					
			// Discrete electronics response sampled in time.
			std::vector<double> k(nKernel, 0.0);    
			for (int ik = 0; ik < nKernel; ++ik) {
					k[ik] = gmct_signalShape2(ik * dt_ns, 1.0);
			}

			// Convolve input waveform with electronics response (signal shaping).
			for (int i = 0; i < nBins; ++i) {
				double accConv = 0.0;
				const int jmax = std::min(i, nKernel - 1);
				for (int j = 0; j <= jmax; ++j) {
						accConv += x[i - j] * k[j];
				}

					// Multiply by dt to approximate continuous convolution integral.
				out[i] = accConv * dt_ns;
			}
		}; // end of apply electronics

		// Left-end and right-end waveform after transmission losses
		// and electronics shaping.
		apply_step10_electronics_vec(wfLeft,  wfElecL);
		apply_step10_electronics_vec(wfRight, wfElecR);

		//============================================================================================//
		//                                add colored noise using FFT/IFFT                            //
		//============================================================================================//
		// Add electronic noise to the shaped waveform using an FFT noise template.
		if (m_enableFFTNoise.value() && m_fftNoiseReady)
		{
			// Use the same samples size as the analog waveform for the noise.
			const int nBins = nWaveformBins;
			// Analog sampling step sets the waveform Nyquist frequency for the noise generation.
			const double dt_ns = waveformBinSize_ns;

			// Generate one colored-noise for the left channel
			// and one for the right channel
			auto noiseL = dchfft::makeFFTNoiseSamples(
				nBins, dt_ns, m_fftMag, m_fftMaxFreq, myRandom, m_noiseRemoveDC.value());
			auto noiseR = dchfft::makeFFTNoiseSamples(
				nBins, dt_ns, m_fftMag, m_fftMaxFreq, myRandom, m_noiseRemoveDC.value());

			// Rescale the generated noise to the requested detector/electronics noise amplitude.
			const double scale = (m_fftAmpNorm != 0.0) ? (m_noiseScale.value() / m_fftAmpNorm) : 0.0;

			// Add the scaled colored noise sample-by-sample to the electronics-shaped waveform.
			for (int i = 0; i < nBins; ++i) {
				wfElecL[i] += scale * noiseL[i];
				wfElecR[i] += scale * noiseR[i];  
			}

			// Fill ONE noise histogram for debugging (only first cell/event)
			if (m_create_debug_histos.value() && !m_waveformFilled) {
				hNoise->Reset();
				for (int ibin = 1; ibin <= nBins; ++ibin) {
					hNoise->SetBinContent(ibin, scale * noiseL[ibin - 1]);
				}
			}
		} // end of adding noise

		//================================================================================================//
		//                     Analog(L/R) -> Digitized(L/R) using the SAME time bin:                     //
		//================================================================================================//
		// Convert analog waveforms to discrete ADC samples
		std::vector<uint16_t> adcL, adcR;
		std::vector<float> digiL_mV, digiR_mV;

		if (m_doWaveformDigitization.value())
		{
			// Maximum ADC code defined by the number of bits
			const int adcMax = (m_adcBits.value() >= 1) ? ((1 << m_adcBits.value()) - 1) : 0;
			// Voltage step per ADC count
			const double lsb_mV = m_adcLSB_mV.value();
			// Baseline offset representing electronics pedestal.
			const double base_mV = m_adcBaseline_mV.value();
			// Electronics polarity defines signal sign convention
			const int pol = (m_adcPolarity.value() >= 0) ? +1 : -1;

			// ADC sampling period defines discrete time grid of digitization.
			const double dtDigi_ns = m_adcSamplePeriod_ns.value();
			// Number of digitized samples within the waveform time window.
			const int nBinsDigi = std::max(
					1, static_cast<int>(std::lround((wfTmax - wfTmin) / dtDigi_ns)));

			// Reserve memory for ADC codes and reconstructed voltages.
			adcL.reserve(nBinsDigi);
			adcR.reserve(nBinsDigi);
			digiL_mV.reserve(nBinsDigi);
			digiR_mV.reserve(nBinsDigi);

			// Loop over ADC sampling times.
			for (int ibin = 0; ibin < nBinsDigi; ++ibin) {
				const double t_ns = wfTmin + (ibin + 0.5) * dtDigi_ns;

				int idxWaveform = static_cast<int>(std::floor((t_ns - wfTmin) / waveformBinSize_ns));
				if (idxWaveform < 0) idxWaveform = 0;
				if (idxWaveform >= nWaveformBins) idxWaveform = nWaveformBins - 1;

				// Voltage at this sampling point
				const double sampledL_V = wfElecL[idxWaveform];
				const double sampledR_V = wfElecR[idxWaveform];

				// Convert voltage to mV, apply polarity and add baseline offset.
				const double vinL_mV = base_mV + pol * (1000.0 * sampledL_V);
				const double vinR_mV = base_mV + pol * (1000.0 * sampledR_V);

				// Quantize input voltage into ADC integer code.
				int codeL = static_cast<int>(vinL_mV / lsb_mV);
				int codeR = static_cast<int>(vinR_mV / lsb_mV);

				// Saturate ADC codes within allowed dynamic range.
				if (codeL < 0) codeL = 0;
				if (codeL > adcMax) codeL = adcMax;
				if (codeR < 0) codeR = 0;
				if (codeR > adcMax) codeR = adcMax;

				// Store digitized ADC counts.
				adcL.push_back(static_cast<uint16_t>(codeL));
				adcR.push_back(static_cast<uint16_t>(codeR));
				digiL_mV.push_back(static_cast<float>(codeL * lsb_mV));
				digiR_mV.push_back(static_cast<float>(codeR * lsb_mV));
			}
		} // end of digitizating alanoge waveform
		//---------------------------------------------------------------------------------//

		// Store one representative waveform per event
		// to monitor the full signal chain without histogram overlap.
		if (isDebugCell)
		{
			if (hWaveform) {
				hWaveform->Reset();
				for (int ibin = 1; ibin <= nWaveformBins; ++ibin) {
					hWaveform->SetBinContent(ibin, wfCenter[ibin - 1]);  
				}  
			}
			if (hWaveformL) {
				hWaveformL->Reset();
				for (int ibin = 1; ibin <= nWaveformBins; ++ibin) {
					hWaveformL->SetBinContent(ibin, wfLeft[ibin - 1]);  
				}  
			}
			if (hWaveformR) {
				hWaveformR->Reset();
				for (int ibin = 1; ibin <= nWaveformBins; ++ibin) {
					hWaveformR->SetBinContent(ibin, wfRight[ibin - 1]);  
				}
			}
			if (hWaveformElecL) {
				hWaveformElecL->Reset();
				for (int ibin = 1; ibin <= nWaveformBins; ++ibin) {
					hWaveformElecL->SetBinContent(ibin, wfElecL[ibin - 1]);
				}
			}
			if (hWaveformElecR) {
				hWaveformElecR->Reset();
				for (int ibin = 1; ibin <= nWaveformBins; ++ibin) {
					hWaveformElecR->SetBinContent(ibin, wfElecR[ibin - 1]);
				}
			}
			if (hWaveformDigiL && !digiL_mV.empty() &&
			(int)digiL_mV.size() == hWaveformDigiL->GetNbinsX()) 
			{
				hWaveformDigiL->Reset();
				for (int ibin = 1; ibin <= hWaveformDigiL->GetNbinsX(); ++ibin) {
					hWaveformDigiL->SetBinContent(ibin, digiL_mV[ibin - 1]);
				}
			}
			if (hWaveformDigiR && !digiR_mV.empty() &&		  
			(int)digiR_mV.size() == hWaveformDigiR->GetNbinsX()) {
				hWaveformDigiR->Reset();
				for (int ibin = 1; ibin <= hWaveformDigiR->GetNbinsX(); ++ibin) {
					hWaveformDigiR->SetBinContent(ibin, digiR_mV[ibin - 1]);
				}
			}
			m_waveformFilled = true;  
		} // End of debugging	  
		//--------------------------------------------------------------------------------------//

		const float t0_ns = static_cast<float>(t_hit_ns);

		// left side digi waveform:
    	if (m_doWaveformDigitization.value() && !adcL.empty()) {
		
			auto DigiwfL = out_Digitized_waveformsL.create();
			DigiwfL.setCellID(kv.first);
			DigiwfL.setTime(t0_ns);
			DigiwfL.setInterval(m_adcSamplePeriod_ns.value());
			DigiwfL.setCharge(m_adcLSB_mV.value());

			for (uint16_t code : adcL) {
				DigiwfL.addToAdcCounts(code);
			}

			//DigiwfL.setHit(oDCHdigihit);
		}
	
		// left side digi waveform:
		if (m_doWaveformDigitization.value() && !adcR.empty()) {

			auto DigiwfR = out_Digitized_waveformsR.create();
			DigiwfR.setCellID(kv.first);
			//DigiwfR.setSide(0);                // convention: 0 = Right
			DigiwfR.setTime(t0_ns);
			DigiwfR.setInterval(m_adcSamplePeriod_ns.value());
			DigiwfR.setCharge(m_adcLSB_mV.value());

			for (uint16_t code : adcR) {
			DigiwfR.addToAdcCounts(code);
			}
		}

	} // End of the Cell loop:

	// return the digitized hits, truth links, and left/right waveform collections for this event
	return std::make_tuple(
		std::move(out_Digitized_waveformsL),
		std::move(out_Digitized_waveformsR)
	);

} // End of the operator Block

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       finalize       //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
StatusCode WaveformFromWireTracker::finalize() 
{
	if (m_polya) {
		delete m_polya;
		m_polya = nullptr;
	}
	return StatusCode::SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       PrintConfiguration       ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
void WaveformFromWireTracker::PrintConfiguration(std::ostream& io) const
{
	io << "DCHdigi will use the following components:\n";
	io << "\tGeometry Service: " << m_geoSvcName.value().c_str() << "\n";
	io << "\tUID Service: " << m_uidSvcName.value().c_str() << "\n";
	io << "\tDetector name: " << m_detectorName.value().c_str() << "\n";
	io << "\t\t|--Volume bitfield: " << m_decoder->fieldDescription().c_str() << "\n";
	io << "\t\t|--Number of layers: " << dch_data->database.size() << "\n";
	io << "\tCreate debug histograms (THistSvc): "
		<< (m_create_debug_histos.value() ? "true" : "false") << "\n";
	return;
}
