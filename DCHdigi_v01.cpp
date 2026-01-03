#include "DCHdigi_v01.h"
#include "BetheBloch.h"
#include "DCHPhysicsTools.h"
#include "XTRELTIME.h"
#include "DCHFFTNoise.h"

#include "TMath.h"
#include "TProfile.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TSpline.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TDirectory.h"
#include "TVirtualFFT.h"
#include <numeric> 
#include "TParameter.h"

// STL
#include <algorithm>
#include <iostream>
#include <sstream>
#include <random>
#include <cmath>
#include <unordered_map>
#include "extension/MutableSenseWireHit.h"

#include "DD4hep/Detector.h"
#include "DD4hep/DetElement.h"
#include "DD4hep/Volumes.h"


namespace {

  // Read a "vector-like" object stored as TGraph or TH1 into std::vector<double>.
  // For your file, ROOT browser shows green icons -> very likely TGraph.
  
  std::vector<double> readVectorFromObj(TObject* obj, bool useX) {
  std::vector<double> out;
  if (!obj) return out;

  if (auto* g = dynamic_cast<TGraph*>(obj)) {
    const int n = g->GetN();
    out.reserve(n);
    for (int i = 0; i < n; ++i) {
      double x = 0.0, y = 0.0;
      g->GetPoint(i, x, y);
      out.push_back(useX ? x : y);
    }
    return out;
  }

  if (auto* h = dynamic_cast<TH1*>(obj)) {
    const int n = h->GetNbinsX();
    out.reserve(n);
    for (int i = 1; i <= n; ++i) {
      out.push_back(useX ? h->GetBinCenter(i) : h->GetBinContent(i));
    }
    return out;
  }

  return out;
}


  // Read a scalar stored as TGraph(1 point) or TH1(1 bin) or TParameter<double>
  double readScalarFromObj(TObject* obj) {
    if (!obj) return 0.0;

    if (auto* p = dynamic_cast<TParameter<double>*>(obj)) {
      return p->GetVal();
    }

    if (auto* g = dynamic_cast<TGraph*>(obj)) {
      double x=0, y=0;
      if (g->GetN() > 0) {
        g->GetPoint(0, x, y);
        return y;
      }
    }

    if (auto* h = dynamic_cast<TH1*>(obj)) {
      if (h->GetNbinsX() >= 1) return h->GetBinContent(1);
    }

    return 0.0;
  }

  inline int float_cmp(double a, double b, double eps=1e-9) {
    if (std::fabs(a-b) < eps) return 0;
    return (a > b) ? 1 : -1;
  }

} // end anon namespace



///////////////////////////////////////////////////////////////////////////////////////
//////////////////////       DCHdigi_v01 constructor       ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
// -- KeyValues("name of the variable that holds the name of the collection exposed in the python steering file",
// {"default name for the collection"}),
DCHdigi_v01::DCHdigi_v01(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc,
                       {
                           KeyValues("DCH_simhits", {""}),
                           KeyValues("HeaderName", {"EventHeader"}),
                       },
                       {KeyValues("DCH_DigiCollection", {"DCH_DigiCollection"}),
                        KeyValues("DCH_DigiSimAssociationCollection", {"DCH_DigiSimAssociationCollection"})}) {
  m_geoSvc = serviceLocator()->service(m_geoSvcName);
  m_uidSvc = serviceLocator()->service(m_uidSvcName);
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       initialize       ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
StatusCode DCHdigi_v01::initialize() {
  if (!m_uidSvc)
    ThrowException("Unable to get UniqueIDGenSvc");

  if (0 > m_z_resolution.value())
    ThrowException("Z resolution input value can not be negative!");

  if (0 > m_xy_resolution.value())
    ThrowException("Radial (XY) resolution input value can not be negative!");

  //-----------------
  // Retrieve the subdetector
  std::string DCH_name(m_DCH_name.value());
  if (0 == m_geoSvc->getDetector()->detectors().count(DCH_name)) {
    ThrowException("Detector <<" + DCH_name + ">> does not exist.");
  }

  dd4hep::DetElement DCH_DE = m_geoSvc->getDetector()->detectors().at(DCH_name);

  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////  retrieve data extension     //////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  this->dch_data = DCH_DE.extension<dd4hep::rec::DCH_info>();
  if (not dch_data->IsValid())
    ThrowException("No valid data extension was found for detector <<" + DCH_name + ">>.");

  ///////////////////////////////////////////////////////////////////////////////////

  //-----------------
  // Retrieve the readout associated with the detector element (subdetector)
  dd4hep::SensitiveDetector dch_sd = m_geoSvc->getDetector()->sensitiveDetector(DCH_name);
  if (not dch_sd.isValid())
    ThrowException("No valid Sensitive Detector was found for detector <<" + DCH_name + ">>.");

  dd4hep::Readout dch_readout = dch_sd.readout();
  // set the cellID decoder
  m_decoder = dch_readout.idSpec().decoder();

  // Read DCH half-length from XML constant (DCH_standalone_o1_v02.xml)
  dd4hep::Detector& det = dd4hep::Detector::getInstance();

  if (!det.constantAsDouble("DCH_half_length_total")) {
  error() << "Missing DD4hep constant: DCH_half_length_total" << endmsg;
  return StatusCode::FAILURE;
  }

  m_halfChamberLength_mm = det.constantAsDouble("DCH_half_length_total"); // returns in mm
  info() << "[DCHdigi] m_halfChamberLength_mm = " << m_halfChamberLength_mm << " mm" << endmsg;

  //Have a look to the cell in the form of a grid:
  gROOT->SetBatch(true);
  const int firstLayer = 80;
  const int nLayers    = 10;
  const int firstPhi   = 0;
  const int nPhi       = 10;

  std::vector<double> xs, ys;
  xs.reserve(nLayers * nPhi);
  ys.reserve(nLayers * nPhi);

  const TVector3 origin(0.0, 0.0, 0.0);
  for (int il = 0; il < nLayers; ++il) {
	  int ilayer = firstLayer + il;
	  for (int ip = 0; ip < nPhi; ++ip) {
		  int nphi = firstPhi + ip;
		  TVector3 v = dch_data->Calculate_hitpos_to_wire_vector(ilayer, nphi, origin);
		  TVector3 p = origin + v;   // point on the wire (closest to origin)
		  xs.push_back(p.X());
		  ys.push_back(p.Y());
	  }
  }
  gSenseWireXY = new TGraph((int)xs.size(), xs.data(), ys.data());
  gSenseWireXY->SetName("grSenseWireXY_10x10");
  gSenseWireXY->SetTitle("Sense-wire positions (10x10 window);x [cm];y [cm]");
  gSenseWireXY->SetMarkerStyle(20);
  gSenseWireXY->SetMarkerSize(0.9);
  
  ///////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////  initialize x-t relation  //////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
 
  m_xt2d = std::make_unique<DCHXT2DLUT>();
  const bool ok = m_xt2d->load(m_xtFileName.value(), "xt_mean", "xt_error", 240, 240);
  if (!ok) {
	  error() << "[XT2D] Failed to load/build 2D XT LUT from " << m_xtFileName.value()
		  << " (need TGraph2D xt_mean and xt_error)" << endmsg;
	  return StatusCode::FAILURE;
  }

  info() << "[XT2D] Built regular-grid XT LUT from " << m_xtFileName.value()
     	  << " rmax~" << m_xt2d->rmaxCm() << " cm"
	  << " grid=" << 240 << "x" << 240 << endmsg;

  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////   Use Noise from par.root  ////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  if (m_enableFFTNoise.value()) {
	  std::string err;
	  const bool ok = dchfft::loadFFTNoiseTemplate(
		    	  m_noiseFileName.value(),
		    	  m_noiseDirName.value(),
		    	  m_fftMag, m_fftAmpNorm, m_fftMaxFreq, m_fftSize,
		    	  &err
			  );

	  if (!ok) {
	      	  error() << "[FFTNoise] " << err << endmsg;
	      	  return StatusCode::FAILURE;
	  }

	  m_fftNoiseReady = true;
	  info() << "[FFTNoise] Loaded fft_noise:"
	 	  << " size=" << m_fftSize
	 	  << " maxFreq=" << m_fftMaxFreq
	 	  << " ampNorm=" << m_fftAmpNorm << endmsg;
  }

  /////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////  initialize Polya gain  ///////////////////////////
  /////////////////////////////////////////////////////////////////////////////////
  // Polya PDF for gas gain:
  // P(Q) ~ ((1+theta)/(G * Gamma(theta+1))) *
  //        [Q*(1+theta)/G]^theta * exp( - Q*(1+theta)/G )
  //
  // Parameters:
  //   [0] normalization (kept at 1; ROOT will renormalize as needed)
  //   [1] mean gas gain G
  //   [2] theta (Polya shape parameter)
  //
  // Range (1e3, 1e7) is taken from the original GMCT implementation
  // and can be tuned later if needed.
  m_polya = new TF1(
      "Polya",
      "((1.+[2])/([1]*TMath::Gamma([2]+1))) * pow(x*(1.+[2])/[1],[2]) * exp(-x*(1.+[2])/[1])",
      1.e3, 1.e7);

  m_polya->SetParameter(0, 1.0);
  m_polya->SetParameter(1, m_GasGain.value());
  m_polya->SetParameter(2, m_PolyaTheta.value());
  //-----------------------------------------------------------------------------------//

  std::stringstream ss;
  PrintConfiguration(ss);
  info() << ss.str().c_str() << endmsg;
  if (m_create_debug_histos.value()) {
    hDpw = new TH1D("hDpw", "Distance from sim-hit to the wire, in cm", 100, 0, 1);
    hDpw->SetDirectory(0);
    hDww = new TH1D(
        "hDww",
        "Distance from digi-hit to the wire, in cm. Should be zero because digi-hit central position lies on the wire",
        100, 0, 1);
    hDww->SetDirectory(0);
    hSz = new TH1D("hSz", "Smearing along the wire, in cm", 100, 0, 5 * m_z_resolution.value());
    hSz->SetDirectory(0);
    hSxy = new TH1D("hSxy", "Smearing perpendicular the wire, in cm", 100, 0, 5 * m_xy_resolution.value());
    hSxy->SetDirectory(0);
  }

//-------------------------------------------New Hist definition--------------------------------------//
  //Histogram for path length:
  hPathLength = new TH1F("hPathLength", "Path Length per Hit;Path Length [mm];Entries", 10, 0, 15);
  hPathLength->SetDirectory(0);
  hPLvsnC = new TH2F("hPLvsnC","Ionized Length vs nCell; PathLength [mm]; nCell", 20, 0, 15,  20, 0, 200);
  hPLvsnC->SetDirectory(0);

  //Histogram for x, y and z positions:
  hX = new TH1F("hX", "Hit X positions;X (mm);Counts", 50, 0, 2000);
  hX->SetDirectory(0);
  hY = new TH1F("hY", "Hit Y positions;Y (mm);Counts", 50, 0, 2000);
  hY->SetDirectory(0);
  hZ = new TH1F("hZ", "Hit Z positions;Z (mm);Counts", 50, 0, 2000);
  hZ->SetDirectory(0);
  hXYZ = new TH3F("hXYZ", "Hit positions;X (mm);Y (mm);Z (mm)", 50, 0, 2000, 50, 0, 2000, 50, 0, 2000);
  hXYZ->SetDirectory(0);

  heDep = new TH1F("eDp", "Energy Deposit", 50, 0, 10);
  heDep->SetDirectory(0);

  hR = new TH1F("hR", "Distance(mm)", 50, 0, 2000);
  hR->SetDirectory(0);
  hRvsZ = new TH2F("hRvsZ","R vs z; z [mm]; R [mm]", 50, 0, 2400,  50, 0, 2400);
  hRvsZ->SetDirectory(0);

  hNcl_perStep = new TH1F("hNcl_perStep", "n_{cl} per sim-step; n_{cl}; entries", 50, 0, 50);
  hNcl_perStep->SetDirectory(0);
  hNe_perStep = new TH1F("hNe_perStep", "N_{e} per sim-step; Cluster size; entries", 10, 0, 10);
  hNe_perStep->SetDirectory(0);
  hClSpacing_mm = new TH1F("hClSpacing_mm", "Inter-cluster spacing; spacing [mm]; entries", 50, 0, 6);
  hClSpacing_mm->SetDirectory(0);

  hTd = new TH1F  ("hTd", "Drift Time; t (ns); Entries", 50, 0, 1000);
  hTd->SetDirectory(0);
  hXT = new TH2F ("hXT", "x-t relation; r (cm); t (ns)", 50, 0, 1.0, 50 , 0, 3500);
  hXT->SetDirectory(0);
  if (m_create_debug_histos.value()) {
  m_grXT2D = new TGraph2D();
  m_grXT2D->SetName("grXT2D");
  m_grXT2D->SetTitle("drift time;X (cm);Y (cm);t (ns)");
  }

  hV = new TH1F("hV", "Drift velocity; V (cm/#mu s); entries", 50, 0, 10);
  hV->SetDirectory(0);
  pVvsE = new TProfile("pVvsE","Electron Drift Velocity vs Electric Field;<E> (V/cm);<v_{d}> (cm/#mu s)",
  50, 0, 5000);
  pVvsE->SetDirectory(0);
  hVvsR = new TH2F ("hVvsR", "Drift Velocity vs Radius; r (cm); V (cm/#mu s)", 50, 0, 1, 50 , 0, 10);
  hVvsR->SetDirectory(0);
  
  hAvalancheQ = new TH1F("hAvalancheQ", "Avalanche per Charge; Q  (arb. units); Entries", 100, 0, 150.e4);
  hAvalancheQ->SetDirectory(0);

  //  Debug: Analog waveform (Point 8)
  const double wfTmin = 0.0;
  const double wfTmax = 2000.0; //m_waveformTimeWindow_ns.value();
  const double dt_ns  = 0.5;
  const int nBinsWave  = static_cast<int>(std::lround((wfTmax - wfTmin) / dt_ns));
  hWaveform = new TH1F("hWaveform", "Analog Waveform; t (ns); Amplitude", nBinsWave, wfTmin, wfTmax);
  hWaveform->SetDirectory(0);
  hWaveformL = new TH1F("hWaveformL", "Waveform at wire end 1; t [ns]; Amplitude", nBinsWave, wfTmin, wfTmax);
  hWaveformL->SetDirectory(0);
  hWaveformR = new TH1F("hWaveformR","Waveform at wire end 2; t [ns]; Amplitude",nBinsWave, wfTmin, wfTmax);
  hWaveformR->SetDirectory(0);
  hWaveformElecL = new TH1F("hWaveformElecL", 
		  "Waveform after electronics (end 1); t [ns]; Amplitude",
		  nBinsWave, wfTmin, wfTmax);
  hWaveformElecL->SetDirectory(0);
  hWaveformElecR = new TH1F("hWaveformElecR",
                          "Waveform after electronics (end 2); t [ns]; Amplitude",
                          nBinsWave, wfTmin, wfTmax);
  hWaveformElecR->SetDirectory(0);

  // Digitized (reconstructed) waveform in mV
  hWaveformDigiL = new TH1F("hWaveformDigiL","Digitized waveform at wire end 1; t [ns]; ADC-reco [mV]",
                          nBinsWave, wfTmin, wfTmax);
  hWaveformDigiL->SetDirectory(0);

  hWaveformDigiR = new TH1F("hWaveformDigiR","Digitized waveform at wire end 2; t [ns]; ADC-reco [mV]",
                          nBinsWave, wfTmin, wfTmax);
  hWaveformDigiR->SetDirectory(0);

  // Optional: raw ADC code
  hADCL = new TH1F("hADCL", "ADC codes end 1; t [ns]; ADC code",
                 nBinsWave, wfTmin, wfTmax);
  hADCL->SetDirectory(0);

  hADCR = new TH1F("hADCR", "ADC codes end 2; t [ns]; ADC code",
                 nBinsWave, wfTmin, wfTmax);
  hADCR->SetDirectory(0);

  hSignalPulse = nullptr;

  m_clusterTree = new TTree("tClusters", "Cluster Count");
  m_clusterTree->SetDirectory(nullptr);
  m_clusterTree->Branch("Ncl_step", &m_cluster_Ncl_step, "Ncl_step/F");

//-----------------------------------------------End-------------------------------------------------//

  return StatusCode::SUCCESS;
}

std::tuple<std::mt19937_64, TRandom3> DCHdigi_v01::CreateRandomEngines(const edm4hep::EventHeaderCollection& headers) const {
  auto engine_seed = m_uidSvc->getUniqueID(headers, this->name());
  auto rng_engine = std::mt19937_64(engine_seed);
  auto random_seed = m_uidSvc->getUniqueID(headers, this->name()+"_1");
  auto myRandom = TRandom3(random_seed);
  // advance internal state to minimize possibility of creating correlations
  rng_engine.discard(10);
  for (int i = 0; i < 10; ++i)
    myRandom.Rndm();
  return {rng_engine, myRandom};
}

///////////////////////////////////////////////////////////////////////////////////////

std::tuple<extension::SenseWireHitCollection, extension::SenseWireHitSimTrackerHitLinkCollection>
DCHdigi_v01::operator()(const edm4hep::SimTrackerHitCollection& input_sim_hits,
                    const edm4hep::EventHeaderCollection&   headers) const {

  /// initialize engines
  auto [rng_engine, myRandom] = this->CreateRandomEngines(headers);
    
  // Gaussian random number generator used for the smearing of the z position, in cm!
  std::normal_distribution<double> gauss_z_cm{0., m_z_resolution.value() * MM_TO_CM};
  // Gaussian random number generator used for the smearing of the xy position, in cm!
  std::normal_distribution<double> gauss_xy_cm{0., m_xy_resolution.value() * MM_TO_CM};

  debug() << "Input Sim Hit collection size: " << input_sim_hits.size() << endmsg;

  // Create the collections we are going to return
  extension::SenseWireHitCollection output_digi_hits;
  extension::SenseWireHitSimTrackerHitLinkCollection output_digi_sim_association;

//==========================================================================================
//                              Changes made by Muhammad saiel                            //
//==========================================================================================
//

  //Varible initailizations:
  edm4hep::Vector3d pos;
  edm4hep::MCParticle mcParticle;
  int pdg;
  float pathLength;
  float R;
  double energyDep;
  int ilayer;
  int nphi;

  // loop over hit collection
  int loop_index = 0;

  for (const auto& input_sim_hit : input_sim_hits) {
    loop_index++;

    pos = input_sim_hit.getPosition();
    energyDep = input_sim_hit.getEDep();
    pathLength = input_sim_hit.getPathLength();
    mcParticle = input_sim_hit.getParticle();
    pdg = mcParticle.getPDG();

    dd4hep::DDSegmentation::CellID cellid = input_sim_hit.getCellID();
    ilayer = this->CalculateLayerFromCellID(cellid);
    nphi = this->CalculateNphiFromCellID(cellid);
    const int hit_ilayer = ilayer;
    const int hit_nphi   = nphi;

    auto hit_position = Convert_EDM4hepVector_to_TVector3(input_sim_hit.getPosition(), MM_TO_CM);
    R = std::sqrt(pos.x*pos.x + pos.y*pos.y);

    //Looking for only MC Particle Hit:
    const bool isPrimaryMuonHit =(std::abs(pdg) == 13) &&
	    !input_sim_hit.isProducedBySecondary();
    bool isDeltaElectronHit =(std::abs(pdg) == 11) &&
	    input_sim_hit.isProducedBySecondary() &&
	    IsParticleCreatedInsideDriftChamber(mcParticle);
    if (!isPrimaryMuonHit) continue;
    //if (isDeltaElectronHit) continue;
    //if(input_sim_hit.getPathLength()<=1.0) continue;
  
    std::cout <<"Hits No: " << loop_index<<"\t"
	    <<"Layer: "<<ilayer<<"\t"
	    <<"Cell Number: "<< nphi <<"\t"
	    << "pathLength: " << input_sim_hit.getPathLength() << " mm"<<"\t"
	    << std::endl;

//=======================================================================================
//                                Calculate N_cluster per cm                           //
//=======================================================================================
    
    const double L_mm = input_sim_hit.getPathLength();
    if (L_mm <=0) continue;
    const double l_cm = 0.1*L_mm;

    const double bg = compute_beta_gamma(mcParticle);
    static const double Ne_mean = mean_cluster_size_He();
    const double lambda = get_dNcldx_per_cm(bg)/Ne_mean;
    double mu = lambda * l_cm;

    //Generate Cluster per cm:
    std::poisson_distribution<int> pois(mu);
    const int Ncl_step = (std::max(0, pois(rng_engine)))/l_cm;
    hNcl_perStep->Fill(Ncl_step);

    //Generate cluster position inside the loop:
    auto cl_positions_mm = generate_cluster_positions_mm(L_mm, lambda, rng_engine);
    for (size_t i = 1; i < cl_positions_mm.size(); ++i)
    {
	    hClSpacing_mm->Fill(cl_positions_mm[i] - cl_positions_mm[i-1]);
    }

    //Generate Cluster size inside the loop:
    int n_clusters_in_step = cl_positions_mm.size();
    int Ne;
    std::vector<int> electrons_per_cluster;
    for (int i = 0; i < n_clusters_in_step; ++i)
    {
	    Ne = sample_cluster_size(rng_engine);
	    electrons_per_cluster.push_back(Ne); 
    }
    for (auto ne : electrons_per_cluster)
    {
	    hNe_perStep->Fill(ne);
    }
    //std::cout<<"No Crash after clustering, you are good to go:"<<std::endl;

    hPathLength->Fill(pathLength);
    hPLvsnC->Fill(pathLength, nphi);
    hX->Fill(pos.x);
    hY->Fill(pos.y);
    hZ->Fill(pos.z);
    hXYZ->Fill(pos.x, pos.y, pos.z);
    heDep->Fill(energyDep);
    hR->Fill(R);
    hRvsZ->Fill(pos.z, R);

    if (m_clusterTree) {
	    m_cluster_Ncl_step = Ncl_step;
	    m_clusterTree->Fill();
    }
//-------------------------------end of the chages----------------------------------//

    // -------------------------------------------------------------------------
    //      calculate hit position projection into the wire
    TVector3 hit_to_wire_vector = this->dch_data->Calculate_hitpos_to_wire_vector(ilayer, nphi, hit_position);
    //double distance_hit_wire = hit_to_wire_vector.Mag();
    //double r_cluster_max = 0.0;
    //static std::vector<double> maxHitDistLayer(2000, 0.0);
    //static std::vector<double> maxCluDistLayer(2000, 0.0);
    //static double globalMaxHitDist_cm = 0.0;
    //static double globalMaxCluDist_cm = 0.0;
    //globalMaxHitDist_cm = std::max(globalMaxHitDist_cm, distance_hit_wire);
    //maxHitDistLayer[ilayer] = std::max(maxHitDistLayer[ilayer], distance_hit_wire);

    TVector3 hit_projection_on_the_wire = hit_position + hit_to_wire_vector;
    if (m_create_debug_histos.value()) {
      double distance_hit_wire = hit_to_wire_vector.Mag();
      hDpw->Fill(distance_hit_wire);
    }
    TVector3 wire_direction_ez = this->dch_data->Calculate_wire_vector_ez(ilayer, nphi);

    // -------------------------------------------------------------------------
    //       smear the position

    //       smear position along the wire
    double smearing_z = gauss_z_cm(rng_engine);
    if (m_create_debug_histos.value())
      hSz->Fill(smearing_z);

    hit_projection_on_the_wire += smearing_z * (wire_direction_ez.Unit());
    if (m_create_debug_histos.value()) {
      // the distance from the hit projection and the wire should be zero
      TVector3 dummy_vector = this->dch_data->Calculate_hitpos_to_wire_vector(ilayer, nphi, hit_projection_on_the_wire);
      hDww->Fill(dummy_vector.Mag());
    }

    //       smear position perpendicular to the wire
    double smearing_xy = gauss_xy_cm(rng_engine);
    if (m_create_debug_histos.value())
      hSxy->Fill(smearing_xy);
    float distanceToWire_real = hit_to_wire_vector.Mag();

    // protect against negative values
    float distanceToWire_smeared = std::max(0.0, distanceToWire_real + smearing_xy);

    std::int32_t type = 0;
    std::int32_t quality = 0;
    float eDepError = 0;
    // length units back to mm
    auto positionSW = Convert_TVector3_to_EDM4hepVector(hit_projection_on_the_wire, 1. / MM_TO_CM);
    // auto  directionSW    = Convert_TVector3_to_EDM4hepVector(wire_direction_ez, 1. / MM_TO_CM);
    float distanceToWire = distanceToWire_smeared / MM_TO_CM;

    // The direction of the sense wires can be calculated as:
    //   RotationZ(WireAzimuthalAngle) * RotationX(stereoangle)
    // One point of the wire is for example the following:
    //   RotationZ(WireAzimuthalAngle) * Position(cell_rave_z0, 0 , 0)
    // variables aredefined below
    auto WireAzimuthalAngle = this->dch_data->Get_cell_phi_angle(ilayer, nphi);
    float WireStereoAngle = 0;
    {
      auto l = this->dch_data->database.at(ilayer);
      // radial middle point of the cell at Z=0
      auto cell_rave_z0 = 0.5 * (l.radius_fdw_z0 + l.radius_fuw_z0);
      // when building the twisted tube, the twist angle is defined as:
      //     cell_twistangle    = l.StereoSign() * DCH_i->twist_angle
      // which forces the stereoangle of the wire to have the oposite sign
      WireStereoAngle = (-1.) * l.StereoSign() * dch_data->stereoangle_z0(cell_rave_z0);
    }

//====================================================================================
//                           Drift Time using x-t realtion                          //
//====================================================================================
    const int total_electrons = std::accumulate(electrons_per_cluster.begin(),
                                            electrons_per_cluster.end(), 0);
    // Vectors to store *per-electron* drift time and radius
    std::vector<double> electron_times_ns;
    electron_times_ns.reserve(total_electrons);	//Just to avoid the reallocation
    std::vector<double> electron_r_cm;
    electron_r_cm.reserve(total_electrons);

    // --- Reconstruct track direction from hit momentum ---
    TVector3 track_dir(input_sim_hit.getMomentum().x,
                   input_sim_hit.getMomentum().y,
                   input_sim_hit.getMomentum().z);
    if (track_dir.Mag() > 0.0) {
	    track_dir = track_dir.Unit();
    }
    else {
	    track_dir.SetXYZ(0.0, 0.0, 0.0); // safety
    }

    TVector3 step_start = hit_position;
    const double L_cm = L_mm*MM_TO_CM;

    // In EDM4hep, SimTrackerHit position is usually the midpoint of the step.
    // So the start point is midpoint - (L/2) * direction.
    if (track_dir.Mag() > 0.0 && L_cm > 0.0) {
	    step_start = hit_position - 0.5 * L_cm * track_dir;
    }

    TVector3 u = wire_direction_ez.Unit();      //unit vector along the wire direction:
    TVector3 e1(1.0, 0.0, 0.0);
    e1 -= (e1.Dot(u)) * u;
    if (e1.Mag() < 1e-9) {
	    e1 = TVector3(0.0, 1.0, 0.0);
	    e1 -= (e1.Dot(u)) * u;
    }
    e1 = e1.Unit();
    TVector3 e2 = u.Cross(e1).Unit();
    
    double s_cm = 0.0;
    int Ne_cluster = 0.0;
    double r_cluster_cm = 0.0;
    double t_ns = 0.0;
    double x_cluster_cm = 0.0;
    double y_cluster_cm = 0.0;

    for (size_t icl = 0; icl < cl_positions_mm.size(); ++icl)
    {
	   s_cm = cl_positions_mm[icl]*MM_TO_CM;		//mm -> cm       
	   Ne_cluster = electrons_per_cluster[icl];

	   TVector3 cl_pos = step_start + s_cm * track_dir;
	   TVector3 cl_to_wire = this->dch_data->Calculate_hitpos_to_wire_vector(hit_ilayer, hit_nphi, cl_pos);
	   // Perpendicular displacement to the wire
	   TVector3 v_perp = cl_to_wire - (cl_to_wire.Dot(u)) * u;

	   // TRUE 2D coordinates in cm (wire-centered transverse plane)
	   x_cluster_cm = v_perp.Dot(e1);
	   y_cluster_cm = v_perp.Dot(e2);

	   // Keep r only for existing debug histos (optional)
	   r_cluster_cm = v_perp.Mag();

	   for (int ie = 0; ie < Ne_cluster; ++ie)
	   {
		   /*const double*/ t_ns = m_xt2d->sampleTimeNs(x_cluster_cm, y_cluster_cm, myRandom);
		   electron_times_ns.push_back(t_ns);
		   electron_r_cm.push_back(r_cluster_cm);

		   hTd->Fill(t_ns);
		   hXT->Fill(r_cluster_cm, t_ns);
	   }
    }
//--------------------------------------------End of Drift Time-----------------------------------------//

    /*
    static int nPrint = 0;
if (++nPrint % 2000 == 0) {
  info() << "[CELL-SIZE-SIMPLE] globalMaxHitDist(cm)=" << globalMaxHitDist_cm
         << " globalMaxClusterDist(cm)=" << globalMaxCluDist_cm
         << endmsg;
}
    
    if (r_cluster_max > distance_hit_wire + 0.5) { // 0.5 cm margin (tunable)
  warning() << "[SUSPECT] hitDist=" << distance_hit_wire
            << " maxClusterDist=" << r_cluster_max
            << " ilayer=" << ilayer << " nphi=" << nphi
            << endmsg;
}

info() << "[CELL-CHECK] hitDist(cm)=" << distance_hit_wire
       << "  maxClusterDist(cm)=" << r_cluster_max
       << "  ratio=" << (distance_hit_wire>0 ? r_cluster_max/distance_hit_wire : -1)
       << endmsg;

    info() << "[CHECK] xt_max_dist_cm=" << xt_max_dist_cm
       << " rmax_seen=" << rmax_seen
       << " out=" << n_out << "/" << n_tot
       << endmsg;

double r_min = 1e9, r_max = -1.0;

for (double r : electron_r_cm) {
  r_min = std::min(r_min, r);
  r_max = std::max(r_max, r);
}
std::cout << "[R] r_min(cm)=" << r_min << "  r_max(cm)=" << r_max << std::endl;
*/

//==========================================================================================
//                               Apply Polya distribution                                 //
//==========================================================================================
    std::vector<double> electron_charges;
    electron_charges.reserve(electron_times_ns.size());
    double Qtot = 0.0;

    for (size_t i = 0; i<electron_times_ns.size(); ++i)
    {
	    double q = avalancheCharge(myRandom);
	    electron_charges.push_back(q);
	    Qtot +=q;

	    if (m_create_debug_histos.value()) hAvalancheQ->Fill(q);
    }

//==========================================================================================
//               Build analog waveform by summing the individual Pulses                   //
//==========================================================================================
    if (m_create_debug_histos.value() && 
		    !m_waveformFilled && 
		    !electron_times_ns.empty()) {
       
      if (!hSignalPulse &&
		      !electron_times_ns.empty() &&
		      electron_times_ns.size() == electron_charges.size()) {
        // Take the earliest electron in this hit
        auto it_min = std::min_element(electron_times_ns.begin(), electron_times_ns.end());
        const size_t idx_min = std::distance(electron_times_ns.begin(), it_min);

        const double t0_exact = *it_min;                 // exact drift time (ns), same reference as waveform
        const double q_exact  = electron_charges[idx_min]; // Polya charge actually used in the waveform

        // Use electronics fall time to choose a plotting window
        const double tau_f = m_pulseFallTime_ns.value();   // ns

        // Time range for plotting the pulse:
        //  - before t0_exact we expect pure baseline (singleElectronPulse returns 0)
        //  - after t0_exact we see the full rise and fall
        double tmin = std::max(0.0, t0_exact - 2.0 * tau_f);
        double tmax = t0_exact + 10.0 * tau_f;  // covers basically the full tail

        const int nBinsPulse = 2000; // fine sampling, totally independent of waveform bins

        hSignalPulse = new TH1F("hSignalPulse",
                                "Exact single-electron pulse; t (ns); amplitude (V)",
                                nBinsPulse, tmin, tmax);
        hSignalPulse->SetDirectory(0);

        // Fill with the REAL pulse used by the digitizer:
        //   V(t) = singleElectronPulse(t, t0_exact, q_exact)
        for (int ibin = 1; ibin <= nBinsPulse; ++ibin) {
          const double t = hSignalPulse->GetBinCenter(ibin);
          const double v = singleElectronPulse(t, t0_exact, q_exact);
          hSignalPulse->SetBinContent(ibin, v);
        }
      }
      // ------------------------------------------------------------------
      // END: single-electron pulse
      // ------------------------------------------------------------------

      // Waveform time axis: 0 .. WaveformTimeWindow_ns
      int nBinsWave = hWaveform->GetNbinsX();
      const double wfTmin  = hWaveform->GetXaxis()->GetXmin();
      const double wfTmax  = hWaveform->GetXaxis()->GetXmax();
      const double wfDelta = (wfTmax - wfTmin) / nBinsWave;

      // Clear any previous content (in case)
      for (int ibin = 1; ibin <= nBinsWave; ++ibin) 
      {
        hWaveform->SetBinContent(ibin, 0.0);
      }

      // Loop over time bins and sum pulses from all electrons
      for (int ibin = 1; ibin <= nBinsWave; ++ibin) 
      {
        const double t = hWaveform->GetBinCenter(ibin);  // waveform time in ns

        double V_t = 0.0;  // waveform value at time t

        for (size_t i = 0; i < electron_times_ns.size(); ++i) {
          const double t0_ns = electron_times_ns[i];   // drift time of electron i (ns)
          const double q     = electron_charges[i];    // Polya amplitude

          // Add contribution of this electron to waveform at time t
          V_t += singleElectronPulse(t, t0_ns, q);
        }

        hWaveform->SetBinContent(ibin, V_t);
      }
// ===========================================================================================
//                   propagate waveform to both ends (delay + attenuation)	            //  
// ===========================================================================================
	const double z_hit_mm  = hit_position.Z() / MM_TO_CM;      // must be in mm here
	const double absz_mm   = std::abs(z_hit_mm);

	const double cosStereo = std::cos(WireStereoAngle);
	const double denom     = (std::abs(cosStereo) > 1e-6) ? cosStereo : 1.0;

	double distR_mm = (m_halfChamberLength_mm - z_hit_mm) / denom;
	double distL_mm = (m_halfChamberLength_mm + z_hit_mm) / denom;
	if (distR_mm < 0.0) distR_mm = 0.0;
	if (distL_mm < 0.0) distL_mm = 0.0;

	const double v_mm_per_ns = m_signalSpeed_mm_per_ns.value();
	const double tpropR_ns = distR_mm / v_mm_per_ns;
	const double tpropL_ns = distL_mm / v_mm_per_ns;

	// Attenuation: exp(-d/lambda)
	const bool doAtt = m_enableWireAttenuation.value() &&
		(m_wireAttenuationLength_mm.value() > 0.0);
	const double lambda_mm = m_wireAttenuationLength_mm.value();

	const double attR = doAtt ? std::exp(-distR_mm / lambda_mm) : 1.0;
	const double attL = doAtt ? std::exp(-distL_mm /  lambda_mm) : 1.0;

	hWaveformL->Reset();
	hWaveformR->Reset();

	for (int ibin = 1; ibin <= hWaveform->GetNbinsX(); ++ibin) {
  	const double t = hWaveform->GetBinCenter(ibin);

  	double VL = 0.0;
  	double VR = 0.0;

  	for (size_t i = 0; i < electron_times_ns.size(); ++i) {
    	const double t0 = electron_times_ns[i];
    	const double q  = electron_charges[i];

    	VR += singleElectronPulse(t, t0 + tpropR_ns, q * attR);
    	VL += singleElectronPulse(t, t0 + tpropL_ns, q * attL);
  	}

 	hWaveformL->SetBinContent(ibin, VL);
  	hWaveformR->SetBinContent(ibin, VR);
	}
// ==========================================================================================//
//                          impedance mismatch (reflection/transmission)                     //
//===========================================================================================//
	auto gmct_signalShape2 = [&](double t_ns, double amp) -> double {
		const double tauFall1 = m_elecTauFall1_ns.value();
		const double tauRise  = m_elecTauRise_ns.value();
		const double mix      = m_elecMixFraction.value();
		const double tauFall2 = m_elecTauFall2_ns.value();

		if (t_ns <= 0.0) return 0.0;
		if (tauFall1 <= 0.0 || tauFall2 <= 0.0 || tauRise <= 0.0) return 0.0;

		double sign = mix * std::exp(-t_ns / tauFall1) / tauFall1;
		sign       += (1.0 - mix) * std::exp(-t_ns / tauFall2) / tauFall2;
		sign       *= (1.0 - std::exp(-t_ns / tauRise));

		 // Normalization factor
		 sign       *= (tauFall1 + tauRise) * (tauRise + tauFall2);
		 sign       /= (tauFall1 * tauFall2 + tauRise * tauFall2
				 + tauFall1 * tauRise * mix - tauRise * mix * tauFall2);
		 return sign * amp;
	};

	auto apply_step10_electronics = [&](TH1F* hIn, TH1F* hOut) {
		if (!hIn || !hOut) return;
		hOut->Reset();
		const int nBins = hIn->GetNbinsX();
		const double dt_ns = hIn->GetXaxis()->GetBinWidth(1);

		// --- mismatch scaling
		double scale = m_frontEndGain.value();
		if (m_enableImpedanceMismatch.value()) {
			const double R = m_matchingRes_Ohm.value();
			const double Z = m_tubeImpedance_Ohm.value();
			const double denom = (R + Z);
			if (std::abs(denom) > 0.0) {
				const double reflect = (R - Z) / denom;
				scale *= (1.0 - reflect);
			}
		}

		// Copy waveform into vector, apply mismatch+gain
		std::vector<double> x(nBins, 0.0);
		for (int ib = 1; ib <= nBins; ++ib) {
			x[ib - 1] = hIn->GetBinContent(ib) * scale;
		}

		// If TF disabled, just copy
		if (!m_enableElectronicsTF.value()) {
			for (int ib = 1; ib <= nBins; ++ib) {
				hOut->SetBinContent(ib, x[ib - 1]);
			}
			return;
		}

		// Discrete convolution
		const int nKernel = std::min(nBins,
				std::max(1, int(std::ceil(m_elecKernelLength_ns.value() / dt_ns)) + 1));
		std::vector<double> k(nKernel, 0.0);
		for (int ik = 0; ik < nKernel; ++ik) {
			k[ik] = gmct_signalShape2(ik * dt_ns, 1.0);
		}

		std::vector<double> y(nBins, 0.0);
		for (int i = 0; i < nBins; ++i) {
			double acc = 0.0;
			const int jmax = std::min(i, nKernel - 1);
			for (int j = 0; j <= jmax; ++j) {
				acc += x[i - j] * k[j];
			}

			// dt factor approximates continuous convolution integral
			y[i] = acc * dt_ns;
		}

		for (int ib = 1; ib <= nBins; ++ib) {
			hOut->SetBinContent(ib, y[ib - 1]);
		}
	};

	apply_step10_electronics(hWaveformL, hWaveformElecL);
	apply_step10_electronics(hWaveformR, hWaveformElecR);

//============================================================================================//
//                                add colored noise using FFT/IFFT                            //
//============================================================================================//
                                
	if (m_enableFFTNoise.value() && m_fftNoiseReady) {
		const int nBins = hWaveformElecL->GetNbinsX();
		const double dt_ns = hWaveformElecL->GetXaxis()->GetBinWidth(1);

		// generate noise samples (same noise model applied independently to each end)
		auto noiseL = dchfft::makeFFTNoiseSamples(nBins, dt_ns, m_fftMag, m_fftMaxFreq,
				myRandom, m_noiseRemoveDC.value());
		auto noiseR = dchfft::makeFFTNoiseSamples(nBins, dt_ns, m_fftMag, m_fftMaxFreq,
				myRandom, m_noiseRemoveDC.value());

		// scale: NoiseGenerator.C does: wf += Re[i] * amplitude / fft_amp
		const double scale = (m_fftAmpNorm != 0.0) ? (m_noiseScale.value() / m_fftAmpNorm) : 0.0;

		for (int ib = 1; ib <= nBins; ++ib) {
			hWaveformElecL->SetBinContent(ib, hWaveformElecL->GetBinContent(ib) + scale * noiseL[ib - 1]);
			hWaveformElecR->SetBinContent(ib, hWaveformElecR->GetBinContent(ib) + scale * noiseR[ib - 1]);
		}
	}

//================================================================================================//
//                     Analog(L/R) -> Digitized(L/R) using the SAME time bin:                     //
//================================================================================================//
	
	std::vector<uint16_t> adcL, adcR;
	std::vector<float>    digiL_mV, digiR_mV;   // optional: reconstructed mV for plotting

	if (m_doWaveformDigitization.value())
	{
	      	const int adcMax = (m_adcBits.value() >= 1) ? ((1 << m_adcBits.value()) - 1) : 0;
	      	const double lsb_mV = m_adcLSB_mV.value();
	      	const double base_mV = m_adcBaseline_mV.value();
		const int pol         = (m_adcPolarity.value() >= 0) ? +1 : -1;

	      	// ---------- Left ----------	      
		{
		    	const int nBins = hWaveformElecL->GetNbinsX();
		    	adcL.reserve(nBins);
		    	digiL_mV.reserve(nBins);

		    	for (int ibin = 1; ibin <= nBins; ++ibin) {
			  	const double analog_mV = 1000.0 * hWaveformElecL->GetBinContent(ibin); // V -> mV
			  	const double vin_mV    = base_mV + pol * analog_mV;                // ADC input voltage

			  	int code = static_cast<int>(vin_mV / lsb_mV);                      // trunc like GMCT
			  	if (code < 0) code = 0;
			  	if (code > adcMax) code = adcMax;

			  	adcL.push_back((uint16_t)code);
			  	digiL_mV.push_back((float)(code * lsb_mV));
		    	}
	      	}
	      	// ---------- Right ----------
	      	{
		    	const int nBins = hWaveformElecR->GetNbinsX();
		    	adcR.reserve(nBins);
		    	digiR_mV.reserve(nBins);

		    	for (int ibin = 1; ibin <= nBins; ++ibin) {
			  	const double analog_mV = 1000.0 * hWaveformElecR->GetBinContent(ibin);
			  	const double vin_mV    = base_mV + pol * analog_mV;

			  	int code = static_cast<int>(vin_mV / lsb_mV);
			  	if (code < 0) code = 0;
			  	if (code > adcMax) code = adcMax;

			  	adcR.push_back((uint16_t)code);
			  	digiR_mV.push_back((float)(code * lsb_mV));
		    	}
	      	}
	}
	// fill the histograms:
	if (m_create_debug_histos.value())
	{
		// Reconstructed digitized waveform (mV)
		if (hWaveformDigiL && !digiL_mV.empty() &&
				(int)digiL_mV.size() == hWaveformDigiL->GetNbinsX()) {
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

		// Optional: raw ADC code plots (very useful debug)
	      	if (hADCL && !adcL.empty() && (int)adcL.size() == hADCL->GetNbinsX()) {
		    	hADCL->Reset();
		    	for (int ibin = 1; ibin <= hADCL->GetNbinsX(); ++ibin) {
			  	hADCL->SetBinContent(ibin, adcL[ibin - 1]);
		    	}
	      	}
	      	if (hADCR && !adcR.empty() && (int)adcR.size() == hADCR->GetNbinsX()) {
		    	hADCR->Reset();
		    	for (int ibin = 1; ibin <= hADCR->GetNbinsX(); ++ibin) {
			  	hADCR->SetBinContent(ibin, adcR[ibin - 1]);
		    	}
	      	}
	}

      // Mark that we have stored one example waveform
      m_waveformFilled = true;
    }
//===============================================================================================//
//                       END OF L0 DIGITIZED WAVEFORM FOR IDEA DCH                               //
//===============================================================================================//                       

    double t_drift_ns =0.0;
    if (!electron_times_ns.empty()) {
	    auto it_min = std::min_element(electron_times_ns.begin(), electron_times_ns.end());
	    t_drift_ns = (*it_min <0.0) ? 0.0 : *it_min;
    }
    const double t_sim_ns = input_sim_hit.getTime();
    const double t_hit_ns = t_sim_ns + t_drift_ns;

    extension::MutableSenseWireHit oDCHdigihit;
    oDCHdigihit.setCellID(input_sim_hit.getCellID());
    oDCHdigihit.setType(type);
    oDCHdigihit.setQuality(quality);
    //oDCHdigihit.setTime(input_sim_hit.getTime());
    oDCHdigihit.setTime(t_hit_ns);
    oDCHdigihit.setEDep(input_sim_hit.getEDep());
    oDCHdigihit.setEDepError(eDepError);
    oDCHdigihit.setPosition(positionSW);
    oDCHdigihit.setPositionAlongWireError(m_z_resolution);
    oDCHdigihit.setWireAzimuthalAngle(WireAzimuthalAngle);
    oDCHdigihit.setWireStereoAngle(WireStereoAngle);
    oDCHdigihit.setDistanceToWire(distanceToWire);
    oDCHdigihit.setDistanceToWireError(m_xy_resolution);
    
    output_digi_hits.push_back(oDCHdigihit);

    extension::MutableSenseWireHitSimTrackerHitLink oDCHsimdigi_association;
    oDCHsimdigi_association.setFrom(oDCHdigihit);
    oDCHsimdigi_association.setTo(input_sim_hit);
    output_digi_sim_association.push_back(oDCHsimdigi_association);

  } // end loop over hit collection
  
/////////////////////////////////////////////////////////////////
  return std::make_tuple<extension::SenseWireHitCollection, extension::SenseWireHitSimTrackerHitLinkCollection>(
      std::move(output_digi_hits), std::move(output_digi_sim_association));
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       finalize       //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

void DCHdigi_v01::Create_outputROOTfile_for_debugHistograms() {
  // save current ROOT directory
  TDirectory* currentDir = gDirectory;

  // save the debug histograms in a file
  // file is saved and closed when going out of scope
  {
    auto filename = m_out_debug_filename.value().c_str();
    std::unique_ptr<TFile> ofile{TFile::Open(filename, "recreate")};
    if (!ofile || ofile->IsZombie()) {
      error() << "Error: Could not open file " << filename << std::endl;
      return;
    }
    ofile->cd();
    hDpw->Write();
    hDww->Write();
    hSxy->Write();
    hSz->Write();
    gSenseWireXY->Write();

    //write the created histograms into the debug.root file
    hPathLength->Write();
    hPLvsnC->Write();
    hX->Write();
    hY->Write();
    hZ->Write();
    hXYZ->Write();
    heDep->Write();
    hR->Write();
    hRvsZ->Write();

    hNcl_perStep->Write();
    hNe_perStep->Write();
    hClSpacing_mm->Write();

    hTd->Write();
    hXT->Write();
    if (m_grXT2D) m_grXT2D->Write();
    hV->Write();
    pVvsE->Write();
    hVvsR->Write();

    hAvalancheQ->Write();
    hSignalPulse->Write();
    hWaveform->Write();
    hWaveformL->Write();
    hWaveformR->Write();
    hWaveformElecL->Write();
    hWaveformElecR->Write();
    hWaveformDigiL->Write();
    hWaveformDigiR->Write();
    hV->Write();
    pVvsE->Write();
    hVvsR->Write();

    hAvalancheQ->Write();
    hSignalPulse->Write();
    hWaveform->Write();
    hWaveformL->Write();
    hWaveformR->Write();
    hWaveformDigiL->Write();
    hWaveformDigiR->Write();
    hADCL->Write();
    hADCR->Write();


    if (m_clusterTree) {
	    m_clusterTree->SetDirectory(ofile.get());
            m_clusterTree->Write();
    }


  }


  // Restore previous ROOT directory
  if (currentDir && (not currentDir->IsDestructed()))
    currentDir->cd();
  return;
}

StatusCode DCHdigi_v01::finalize() {
  if (m_create_debug_histos.value()) {
    this->Create_outputROOTfile_for_debugHistograms();
  }

  if (m_xtHelper) {
	  delete m_xtHelper;
	  m_xtHelper = nullptr;
  }
  if (m_polya) {
	  delete m_polya;
	  m_polya = nullptr;
  }


  return StatusCode::SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       ThrowException       ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
void DCHdigi_v01::ThrowException(std::string s) const {
  error() << s.c_str() << endmsg;
  throw std::runtime_error(s);
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       PrintConfiguration       ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
void DCHdigi_v01::PrintConfiguration(std::ostream& io) {
  io << "DCHdigi will use the following components:\n";
  io << "\tGeometry Service: " << m_geoSvcName.value().c_str() << "\n";
  io << "\tUID Service: " << m_uidSvcName.value().c_str() << "\n";
  io << "\tDetector name: " << m_DCH_name.value().c_str() << "\n";
  io << "\t\t|--Volume bitfield: " << m_decoder->fieldDescription().c_str() << "\n";
  io << "\t\t|--Number of layers: " << dch_data->database.size() << "\n";
  io << "\tResolution along the wire (mm): " << m_z_resolution.value() << "\n";
  io << "\tResolution perp. to the wire (mm): " << m_xy_resolution.value() << "\n";
  io << "\tCreate debug histograms: " << (m_create_debug_histos.value() ? "true" : "false") << "\n";
  if (true == m_create_debug_histos.value())
    io << "\t\t|--Name of output file with debug histograms: " << m_out_debug_filename.value() << "\n";

  return;
}

bool DCHdigi_v01::IsParticleCreatedInsideDriftChamber(const edm4hep::MCParticle& thisParticle) const {
  auto vertex = thisParticle.getVertex(); // in mm
  auto vertexRsquared = vertex[0] * vertex[0] + vertex[1] * vertex[1];
  auto vertexZabs = std::fabs(vertex[2]);
  float DCH_halflengh = dch_data->Lhalf / dd4hep::mm;                // 2000;
  float DCH_rin_squared = std::pow(dch_data->rin / dd4hep::mm, 2);   // 350 * 350;
  float DCH_rout_squared = std::pow(dch_data->rout / dd4hep::mm, 2); // 2000 * 2000;
  return (vertexZabs < DCH_halflengh) && (vertexRsquared > DCH_rin_squared) && (vertexRsquared < DCH_rout_squared);
}
