#include "DCHdigi_v01.h"
#include "BetheBloch.h"
#include "DCHPhysicsTools.h"
#include "XTRELTIME.h"

#include "TMath.h"

// STL
#include <algorithm>
#include <iostream>
#include <sstream>
#include <random>
#include <cmath>
#include <unordered_map>
#include "extension/MutableSenseWireHit.h"

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
  
  ///////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////  initialize x-t relation  //////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  info() << "Loading x-t relation from file: " << m_xtFileName.value() <<endmsg;
  m_xtHelper = new XTRELTIME(m_xtFileName.value().c_str());
 
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
  //-------------------------------------Changes made bu Muhammad Saiel---------------------------------//
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

  hMC_perCell = new TH1F("hMC_perCell", "Mean Number of Cluster; Mean_{Cluster}; Entries", 30, 0, 50);
  hMC_perCell->SetDirectory(0);

  hNcl_perStep = new TH1F("hNcl_perStep", "n_{cl} per sim-step; n_{cl}; entries", 50, 0, 50);
  hNcl_perStep->SetDirectory(0);

  hNe_perStep = new TH1F("hNe_perStep", "N_{e} per sim-step; Cluster size; entries", 10, 0, 10);
  hNe_perStep->SetDirectory(0);

  hClSpacing_mm = new TH1F("hClSpacing_mm", "Inter-cluster spacing; spacing [mm]; entries", 50, 0, 6);
  hClSpacing_mm->SetDirectory(0);

  hXT = new TH2F ("hXT", "x-t relation; r (cm); t (ns)", 50, 0, 1, 50 , 0, 400);
  hXT->SetDirectory(0);
  hT = new TH1F ("hT", "Drift Time; t (ns); Entries", 30, 0, 300);
  hT->SetDirectory(0);
  hV = new TH1F ("hV", "Effective Drift Velocity; V (cm/#mu s); Entries", 50, 0, 30);
  hV->SetDirectory(0);
  hVvsR = new TH2F ("hVvsR", "Drift Velocity vs Radius; r (cm); V (cm/#mu s)", 50, 0, 1, 50 , 0, 100);
  hVvsR->SetDirectory(0);
  hDriftTime = new TH1F ("hDriftTime", "Drift Time (t_{sim} - t_{hit}); t_{drift} (ns); Entries", 30, 0, 300);
  hDriftTime->SetDirectory(0);
  hDriftTimeVsLayer = new TH2F ("hDriftTimeVsLayer", "Drift Time VS Layer; Layer; t_{drift} (ns)", 20, 0, 112, 20, 0, 150);
  hDriftTimeVsLayer->SetDirectory(0);

  hAvalancheQ = new TH1F("hAvalancheQ", "Avalanche per Charge; Q  (arb. units); Entries", 100, 0, 150.e3);
  hAvalancheQ->SetDirectory(0);

  //  Debug: Analog waveform (Point 8)
  int nBinsWave = static_cast<int>(m_waveformTimeWindow_ns.value() / m_waveformBinSize_ns.value());
  //if (nBinsWave < 10) nBinsWave = 10;  // safety
  const double wfTmin = 100.0;
  const double wfTmax = 200.0; //m_waveformTimeWindow_ns.value();
  hWaveform = new TH1F("hWaveform", "Analog Waveform; t (ns); Amplitude", nBinsWave, wfTmin, wfTmax);
  hWaveform->SetDirectory(0);

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
 
//------------------------------Changes made by Muhammad saiel----------------------------//

  //Varible initailizations:
  edm4hep::Vector3d pos;
  edm4hep::MCParticle mcParticle;
  int pdg;
  float pathLength;
  float R;
  double energyDep;
  int ilayer;
  int nphi;

  std::unordered_map<dd4hep::DDSegmentation::CellID, CellAcc> acc;
  // loop over hit collection
  int loop_index = 0;
  for (const auto& input_sim_hit : input_sim_hits) {
    loop_index++;

    // This block code is looking for pathlenght of the and postion of the hit
    pos = input_sim_hit.getPosition();
    energyDep = input_sim_hit.getEDep();
    mcParticle = input_sim_hit.getParticle();
    pdg = mcParticle.getPDG();

    const bool isPrimaryMuonHit = 
	    (std::abs(pdg) == 13) &&
	    !input_sim_hit.isProducedBySecondary();
    bool isDeltaElectronHit = 
	    (std::abs(pdg) == 11) &&
	    input_sim_hit.isProducedBySecondary() &&
	    IsParticleCreatedInsideDriftChamber(mcParticle);

    // Apply filter for only MCParticle hits:
    if (!isPrimaryMuonHit) continue;
    if (isDeltaElectronHit) continue;
    //if(input_sim_hit.getPathLength()<=1.0) continue;

    dd4hep::DDSegmentation::CellID cellid = input_sim_hit.getCellID();
    ilayer = this->CalculateLayerFromCellID(cellid);
    nphi = this->CalculateNphiFromCellID(cellid);
    auto hit_position = Convert_EDM4hepVector_to_TVector3(input_sim_hit.getPosition(), MM_TO_CM);

    pathLength = input_sim_hit.getPathLength();
    R = std::sqrt(pos.x*pos.x + pos.y*pos.y);
    std::cout <<"Primary Ionization: " << loop_index<<"\t"
	    <<"Layer: "<<ilayer<<"\t"
	    <<"Cell Number: "<< nphi <<"\t"
	    << "pathLength: " << input_sim_hit.getPathLength() << " mm"<<"\t"
	    // "Position (x,y,z): (" << pos.x << ", " << pos.y << ", " << pos.z<<" ) mm"<<"\t"
	    // "Tracke: R = "<< R <<"\t"
            //<< "PDG: " << mcParticle.getPDG() << "\t"
            //<< "MCParticle ID: " << mcParticle.id()
	    //<< "Hit Energy: " << input_sim_hit.getEDep()<<"GeV"
	    //<<"Path Length/Cell : "<<a.length_mm<<" mm"
	    << std::endl;

    //---------------------Calculate N_cluster per step-------------------------//
    const double L_mm = input_sim_hit.getPathLength();
    if (L_mm <=0) continue;
    const double l_cm = 0.1*L_mm;

    const double bg = compute_beta_gamma(mcParticle);
    const double lambda_per_cm = get_dNcldx_per_cm(bg);
    double mu = lambda_per_cm * l_cm;

    auto& a = acc[cellid];
    a.length_mm += L_mm;
    a.bgL_sum += bg * L_mm;
    a.mu_sum += lambda_per_cm * l_cm;

    std::poisson_distribution<int> pois(mu);
    const int Ncl_step = (std::max(0, pois(rng_engine)))/l_cm;
    hNcl_perStep->Fill(Ncl_step);

    //Generate cluster position inside the loop:
    auto cl_positions_mm = generate_cluster_positions_mm(L_mm, lambda_per_cm, rng_engine);
    int n_clusters_in_step = cl_positions_mm.size();
    for (size_t i = 1; i < cl_positions_mm.size(); ++i)
    {
      hClSpacing_mm->Fill(cl_positions_mm[i] - cl_positions_mm[i-1]);
    }
    
    //Generate Cluster size inside the loop:
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
//--------------------------------------------------------------------------------------------//
//-----------------------------------------end of the chages----------------------------------//
//--------------------------------------------------------------------------------------------//



    // -------------------------------------------------------------------------
    //      calculate hit position projection into the wire
    TVector3 hit_to_wire_vector = this->dch_data->Calculate_hitpos_to_wire_vector(ilayer, nphi, hit_position);
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

    //------------------------------Drift Time------------------------------//
    
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
    if (track_dir.Mag() > 0.0) 
    {
	    track_dir = track_dir.Unit();
    }
    else 
    {
	    track_dir.SetXYZ(0.0, 0.0, 0.0); // safety
    }

    TVector3 step_start = hit_position;
    const double L_cm = 0.1*L_mm;

    // In EDM4hep, SimTrackerHit position is usually the midpoint of the step.
    // So the start point is midpoint - (L/2) * direction.
    if (track_dir.Mag() > 0.0 && L_cm > 0.0) 
    {
	    step_start = hit_position - 0.5 * L_cm * track_dir;
    }

    for (size_t icl = 0; icl < cl_positions_mm.size(); ++icl)
    {
	   const double s_cm       = cl_positions_mm[icl]*0.1; //mm -> cm       
	   const int    Ne_cluster = electrons_per_cluster[icl];
	   TVector3 cl_pos = step_start + s_cm * track_dir;
	   
	   // Distance from this cluster to the sense wire
	   TVector3 cl_to_wire = this->dch_data->Calculate_hitpos_to_wire_vector(ilayer, nphi, cl_pos);
	   const double r_cluster_cm = cl_to_wire.Mag();

	   // For each electron in this cluster, sample a drift time
	   for (int ie = 0; ie < Ne_cluster; ++ie)
	   {
		   double t_ns = electronDriftTime(r_cluster_cm, myRandom);
		   electron_times_ns.push_back(t_ns);
	           electron_r_cm.push_back(r_cluster_cm);

		   if (m_create_debug_histos.value())
		   {
			   hXT->Fill(r_cluster_cm, t_ns);
		   }
	   }
    }

    //------------------------Apply Polya distribution---------------------//
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

    //-------Build analog waveform by summing the individual waveform-----//
    if (m_create_debug_histos.value() && 
		    !m_waveformFilled && 
		    !electron_times_ns.empty()) {
       
      // ------------------------------------------------------------------
      //                        single-electron pulse
      // ------------------------------------------------------------------
      if (!hSignalPulse && !electron_charges.empty()) {
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

      // Mark that we have stored one example waveform
      m_waveformFilled = true;
    }
    // ------------------------------------------------------------------



    double t_drift_ns =0.0;
    if (!electron_times_ns.empty()) {
	    auto it_min = std::min_element(electron_times_ns.begin(), electron_times_ns.end());
	    t_drift_ns = (*it_min <0.0) ? 0.0 : *it_min;
    }

    const double t_sim_ns = input_sim_hit.getTime();
    const double t_hit_ns = t_sim_ns + t_drift_ns;
    
    if (m_create_debug_histos.value()) {
	    for (size_t i = 0; i < electron_times_ns.size(); ++i) {
		    const double t_ns = electron_times_ns[i];
		    const double r_cm = electron_r_cm[i];
		    if (t_ns > 0.0) {
			    const double V = (r_cm / t_ns) * 1000.0; // cm/us
			    hV->Fill(V);
			    hVvsR->Fill(r_cm, V);
		    }
	    }
    }


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

    hMC_perCell->Write();
    hNcl_perStep->Write();
    hNe_perStep->Write();
    hClSpacing_mm->Write();

    hXT->Write();
    hT->Write();
    hV->Write();
    hVvsR->Write();
    hDriftTime->Write();
    hDriftTimeVsLayer->Write();

    hAvalancheQ->Write();
    hSignalPulse->Write();
    hWaveform->Write();

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

