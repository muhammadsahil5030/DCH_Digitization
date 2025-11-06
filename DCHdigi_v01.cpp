#include "DCHdigi_v01.h"
#include "BetheBloch.h"
#include "TMath.h"
#include "TH1F.h"
#include "TFile.h"

// STL
#include <iostream>
#include <sstream>

#include "extension/MutableSenseWireHit.h"

// random number generator (RNG)
#include <random>
#include <unordered_map>
#include <cmath>

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
  //////////////////  initialize Walaa's code for Cluster counting  /////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  debug() << Form("Opening %s...", m_fileDataAlg.value().c_str()) << endmsg;
  if (not IsFileGood(m_fileDataAlg.value()))
    ThrowException("File <<" + m_fileDataAlg.value() + ">> not found.");
  flData = new AlgData();
  flData->Load_file(m_fileDataAlg.value().c_str());
  flData->Load_interp();
  debug() << Form("Opening %s... Done", m_fileDataAlg.value().c_str()) << endmsg;
  ///////////////////////////////////////////////////////////////////////////////////

  std::stringstream ss;
  PrintConfiguration(ss);
  info() << ss.str().c_str() << endmsg;
  if (m_create_debug_histos.value())
  {
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

  hTotPathCell = new TH1F("hTotPathCell", "total Path length/Cell", 50, 0, 200);
  hTotPathCell->SetDirectory(0);

  hBetaGammaCell = new TH1F("hBetaGammaCell", "Value of BetaGamma/Cell", 50, 0, 200);
  hBetaGammaCell->SetDirectory(0);
  
  hMC_perCell = new TH1F("hMC_perCell", "Mean Number of Cluster; Mean_{Cluster}; Entries", 30, 0, 50);
  hMC_perCell->SetDirectory(0);

  hNcl_perStep = new TH1F("hNcl_perStep", "n_{cl} per sim-step; n_{cl}; entries", 20, 0, 50);
  hNcl_perStep->SetDirectory(0);
  hNcl_perStep10 = new TH1F("hNcl_perStep10", "n_{cl} per sim-step; n_{cl}; entries", 20, 0, 50);
  hNcl_perStep10->SetDirectory(0);
  hNcl_perStep2_5 = new TH1F("hNcl_perStep2.5", "n_{cl} per sim-step; n_{cl}; entries", 20, 0, 50);
  hNcl_perStep2_5->SetDirectory(0);
  hNcl_perStep6 = new TH1F("hNcl_perStep6", "n_{cl} per sim-step; n_{cl}; entries", 20, 0, 50);
  hNcl_perStep6->SetDirectory(0);
  hNcl_perCell = new TH1F("hNcl_perCell", "n_{cl} per Cell; N_{cl}; entries", 20, 0, 40);
  hNcl_perCell->SetDirectory(0);

  hNe_perStep = new TH1F("hNe_perStep", "N_{e} per sim-step; Cluster size; entries", 10, 0, 10);
  hNe_perStep->SetDirectory(0);
  hNe_perCell = new TH1F("hNe_perCell", "N_{e} per Cell; Cluster size; entries", 10, 0, 10);
  hNe_perCell->SetDirectory(0);

  hClSpacing_mm = new TH1F("hClSpacing_mm", "Inter-cluster spacing; spacing [mm]; entries", 50, 0, 6);
  hClSpacing_mm->SetDirectory(0);
  hClSpacing_cell = new TH1F("hClSpacing_cell", "Inter-cluster spacing; spacing [mm]; entries", 50, 0, 6);
  hClSpacing_cell->SetDirectory(0);

  hNcl_vs_l = new TH2F("hNcl_vs_l", "n_{cl} vs path length; l [mm]; n_{cl}", 60, 0, 30, 20, 0, 50);
  hNcl_vs_l->SetDirectory(0);
  
  //These histogram show the walaa algorithm Number of clusters
  hNcl_perStep_Walaa = new TH1F("hNcl_perStep_Walaa","N_{cl} per step (Walaa)", 20, 0, 40);
  
  hClSpacing_Walaa_mm = new TH1F("hClSpacing_Walaa_mm","Cluster spacing (Walaa);mm;entries",100,0,2);
  
  hNcl_vs_l_Walaa = new TH2F("hNcl_vs_l_Walaa","N_{cl} vs l (Walaa);l [cm];N_{cl}",50,0,5,50,0,50);

  hNcl_perStep_Walaa->SetDirectory(0);
  hClSpacing_Walaa_mm->SetDirectory(0);
  hNcl_vs_l_Walaa->SetDirectory(0);

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

//------------------------------------------Changes made by M. Saiel---------------------------------//
//Accumulator:
struct CellAcc 
{
  double length_mm = 0.0;
  double bgL_sum = 0.0;
};

// beta*gamma = |p|/m from MC truth
inline double compute_beta_gamma(const edm4hep::MCParticle& mc)
{
  const auto& p = mc.getMomentum();        // GeV
  const double p_mag = std::sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
  const double m     = mc.getMass();       // GeV
  if (m <= 0.0) return 0.0;
  return p_mag / m;
}

// This function is used to calculate the dNdx using BB.h file
double DCHdigi_v01::get_dNcldx_per_cm(double betagamma) const
{
  // 1) reconstruct momentum from βγ and the configured reference mass
  const double m_GeV = m_MassForBB_GeV.value();
  const double p_GeV = betagamma * m_GeV;

  // 2) call your header: bethe_bloch(xp, par)
  double xp[1]  = { p_GeV };                           // GeV (header multiplies by 1e3 internally)
  double par[2] = { m_GeV * 1.0e3,                     // mass in MeV
                    m_MeanExcEnergy_eV.value() * 1e-6  // I in MeV
                  };
  const double dEdx_MeVcm2_per_g = BB::bethe_bloch(xp, par);

  // 3) convert to MeV/cm using A (g/cm^3)
  const double dEdx_MeV_per_cm = dEdx_MeVcm2_per_g * m_GasDensity_g_cm3.value();

  // 4) convert stopping power to *cluster rate* λ [clusters/cm] using W_eff (per cluster)
  const double lambda_per_cm = (dEdx_MeV_per_cm * 1.0e6) / m_W_eff_eV.value();

  // Safety: never negative
  return (lambda_per_cm > 0.0) ? lambda_per_cm : 0.0;
}

// Grouping all function for calculating number of clusters:
namespace
{
// Generate exponentially-distributed cluster positions along a step of length l_mm.
// Output: vector of positions in mm, measured from the start of the step.
inline std::vector<double>
generate_cluster_positions_mm(double l_mm, double lambda_per_cm, std::mt19937_64& rng)
{
  std::vector<double> pos_mm;
  if (l_mm <= 0.0 || lambda_per_cm <= 0.0) return pos_mm;

  const double l_cm = 0.1 * l_mm;
  // Mean spacing = 1/lambda ; parameter of exponential = 1/mean
  std::exponential_distribution<double> expo(lambda_per_cm);

  double s_cm = 0.0;
  while (true) {
    s_cm += expo(rng);           // draw next spacing
    if (s_cm >= l_cm) break;     // stop if we exceed the step length
    pos_mm.push_back(10.0 * s_cm); // store position in mm
  }
  return pos_mm;

}
// --- Experimental cluster size probabilities w(n) for He–isobutane (90/10) ---
// Based on Fischle et al., NIM A301 (1991)
inline const std::vector<double>& w_cluster_He_iC4H10()
{
  static const std::vector<double> w = {
    0.78, 0.12, 0.034, 0.016, 0.0095, 0.006, 0.0044, 0.0034,
    0.0027, 0.0021, 0.0017, 0.0013, 0.0010, 0.0008, 0.0006
  };
  static bool normalized = false;
  if (!normalized) {
    double sum = std::accumulate(w.begin(), w.end(), 0.0);
    for (auto& x : const_cast<std::vector<double>&>(w)) x /= sum;
    normalized = true;
  }
  return w;
}

// --- Draw a cluster size n_e according to experimental probabilities ---
inline int sample_cluster_size(std::mt19937_64& rng)
{
  const auto& w = w_cluster_He_iC4H10();
  std::discrete_distribution<int> dist(w.begin(), w.end());
  int n = dist(rng) + 1; // +1 because bins start at 1
  return n;
}

} //end namespace
//-----------------------------------------End of the changes---------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       operator()       ////////////////////////////////////////
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

//------------------------------------------------------------------------------------------//
//---------------------------Changes by Muhammad Saiel--------------------------------------//
//------------------------------------------------------------------------------------------//
  
  std::map<dd4hep::DDSegmentation::CellID, CellAcc> acc;
  std::cout << "input_sim_hits.size(): " << input_sim_hits.size() << std::endl;
  fflush(stdout);

  int num_layers = this->dch_data->nlayers;
  int layers_per_superlayer = 8;
  int total_cells = 0;
  int total_layers =0;
  int total_superlayers = (num_layers + layers_per_superlayer - 1) / layers_per_superlayer;

  int loop_index = 0;
  for (int ilayer = 1; ilayer <= num_layers; ++ilayer)
  {

    // Check layer exists in database
    if (this->dch_data->database.find(ilayer) == this->dch_data->database.end())
    {
      std::cout << "Warning: ilayer " << ilayer << " not found in dch_data->database." << std::endl;
      continue;
    }

    int ncell = this->dch_data->Get_ncells(ilayer);
    int superlayer = ilayer / layers_per_superlayer + 1;

    total_layers ++;
    total_cells += ncell;

    std::cout << "Superlayer: " << superlayer
              << " | Layer (global): " << ilayer
              << " | Cells in layer: " << ncell << std::endl;

  }// end of the loop on layers:

  std::cout << ">>>>>>>> DCHdigi.cpp is imported >>>>>>>>>";
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "IDEA Drift Chamber Layers and Cells Summary" << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "----------------------------------------------" << std::endl;
  std::cout << "Total number of cells in DCH: " << total_cells << std::endl;
  std::cout << "Total number of layers in DCH: " << total_layers << std::endl;
  std::cout << "Total number of superlayers in DCH: " << total_superlayers << std::endl;
  std::cout << "----------------------------------------------" << std::endl;

  //--------------------------------------End of my codes---------------------------------------//
  //--------------------------------------------------------------------------------------------//

  edm4hep::Vector3d pos;
  edm4hep::MCParticle mcParticle;
  int pdg;
  float pathLength;
  float R;
  double energyDep;
  int ilayer;
  int nphi;

  // loop over hit collection
  for (const auto& input_sim_hit : input_sim_hits) {
    loop_index ++;

    //------------------------------Changes made by Muhammad saiel----------------------------//
    
    // This block code is looking for pathlenght of the and postion of the hit
    pos = input_sim_hit.getPosition();
    energyDep = input_sim_hit.getEDep();
    mcParticle = input_sim_hit.getParticle();
    pdg = mcParticle.getPDG();

    // apply condition for only mouns hits
    if (pdg != -13) continue;
    dd4hep::DDSegmentation::CellID cellid = input_sim_hit.getCellID();
    ilayer = this->CalculateLayerFromCellID(cellid);
    nphi = this->CalculateNphiFromCellID(cellid);
    auto hit_position = Convert_EDM4hepVector_to_TVector3(input_sim_hit.getPosition(), MM_TO_CM);
    pathLength = input_sim_hit.getPathLength();
    R = std::sqrt(pos.x*pos.x + pos.y*pos.y);
    std::cout <<"Index: " << loop_index<<"\t"
	    //<<"Layer: "<<ilayer<<"\t"
	    //<<"Cell Number: "<< nphi <<"\t"
	    //<< "pathLength: " << input_sim_hit.getPathLength() << " mm"<<"\t"
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

    std::poisson_distribution<int> pois(mu);
    const int Ncl_step = std::max(0, pois(rng_engine));

    auto& a = acc[cellid];
    a.length_mm += L_mm;
    a.bgL_sum += bg * L_mm;

    if (m_create_debug_histos.value())
    {
            hNcl_perStep->Fill(Ncl_step);
            if (L_mm > 9.5 && L_mm < 10.5)
            {
                    hNcl_perStep10->Fill(Ncl_step);
            }
            if (L_mm > 2.3 && L_mm < 3.3)
            {
                    hNcl_perStep2_5->Fill(Ncl_step);
            }
            if (L_mm > 5.5 && L_mm < 6.5)
            {
                    hNcl_perStep6->Fill(Ncl_step);
            }
	}
    
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

    //--------------------------------------------------------------------------/

    hPathLength->Fill(pathLength);
    hPLvsnC->Fill(pathLength, nphi);

    hX->Fill(pos.x);
    hY->Fill(pos.y);
    hZ->Fill(pos.z);
    hXYZ->Fill(pos.x, pos.y, pos.z);

    heDep->Fill(energyDep);
    hR->Fill(R);
    hRvsZ->Fill(pos.z, R);
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

    extension::MutableSenseWireHit oDCHdigihit;
    oDCHdigihit.setCellID(input_sim_hit.getCellID());
    oDCHdigihit.setType(type);
    oDCHdigihit.setQuality(quality);
    oDCHdigihit.setTime(input_sim_hit.getTime());
    oDCHdigihit.setEDep(input_sim_hit.getEDep());
    oDCHdigihit.setEDepError(eDepError);
    oDCHdigihit.setPosition(positionSW);
    oDCHdigihit.setPositionAlongWireError(m_z_resolution);
    oDCHdigihit.setWireAzimuthalAngle(WireAzimuthalAngle);
    oDCHdigihit.setWireStereoAngle(WireStereoAngle);
    oDCHdigihit.setDistanceToWire(distanceToWire);
    oDCHdigihit.setDistanceToWireError(m_xy_resolution);
    // For the sake of speed, let the dNdx calculation be optional
    if (m_calculate_dndx.value()) {
      auto [nCluster, nElectrons_v] = CalculateClusters(input_sim_hit, myRandom);
      // to return the total number of electrons within the step, do the following:
      //   int nElectronsTotal = std::accumulate( nElectrons_v.begin(), nElectrons_v.end(), 0);
      //   oDCHdigihit.setNElectronsTotal(nElectronsTotal);
      // to copy the vector of each cluster size to the EDM4hep data extension, do the following:
      for (auto ne : nElectrons_v)
        oDCHdigihit.addToNElectrons(ne);
    }

    output_digi_hits.push_back(oDCHdigihit);

    extension::MutableSenseWireHitSimTrackerHitLink oDCHsimdigi_association;
    oDCHsimdigi_association.setFrom(oDCHdigihit);
    oDCHsimdigi_association.setTo(input_sim_hit);
    output_digi_sim_association.push_back(oDCHsimdigi_association);

  } // end loop over hit collect

//--------------------------------Changes made by Muhammad saiel-------------------------------//

for (const auto& kValue : acc)
{
  const auto& a = kValue.second;
  if (a.length_mm <= 0.0) continue;

  const double bg_cell = a.bgL_sum / a.length_mm;
  const double lambda_cell_per_cm = get_dNcldx_per_cm(bg_cell);
  const double mu_cell = lambda_cell_per_cm * (0.1*a.length_mm);
  if (a.length_mm > 10 && a.length_mm <20)
  {
        hMC_perCell->Fill(mu_cell);
  }

  std::poisson_distribution<int> pois_cell(mu_cell);
  const int Ncl_cell = std::max(0, pois_cell(rng_engine));
  if (m_create_debug_histos.value())
  {
          if (a.length_mm > 10 && a.length_mm <20)
          {
                hNcl_perCell->Fill(Ncl_cell);
          }
  }
//Generate Cluster position per cell:
  auto cl_positions_cell = generate_cluster_positions_mm(a.length_mm, lambda_cell_per_cm, rng_engine);
  int n_clusters_in_cell = cl_positions_cell.size();
  for (size_t cell = 1; cell < cl_positions_cell.size(); ++cell)
    {
      hClSpacing_cell->Fill(cl_positions_cell[cell] - cl_positions_cell[cell-1]);
    }

  //Generate Cluster size per Cell:
  int ne;
  std::vector<int> electrons_per_cluster_cell;
  for (int j = 0; j < n_clusters_in_cell; ++j)
  {
          ne = sample_cluster_size(rng_engine);
          electrons_per_cluster_cell.push_back(ne);

  }
  for (auto n : electrons_per_cluster_cell)
  {
          hNe_perCell->Fill(n);
  }
}// end of the loop for the cluster information per cell:

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

    hTotPathCell->Write();
    hBetaGammaCell->Write();
    
    hMC_perCell->Write();

    hNcl_perStep->Write();
    hNcl_perStep10->Write();
    hNcl_perStep2_5->Write();
    hNcl_perStep6->Write();
    hNcl_perCell->Write();

    hNe_perStep->Write();
    hNe_perCell->Write();

    hClSpacing_mm->Write();
    hClSpacing_cell->Write();

    hNcl_vs_l->Write();

    hNcl_perStep_Walaa->Write();
    hClSpacing_Walaa_mm->Write();
    hNcl_vs_l_Walaa->Write();
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
  io << "\tCluster distributions taken from: " << m_fileDataAlg.value().c_str() << "\n";
  io << "\tResolution along the wire (mm): " << m_z_resolution.value() << "\n";
  io << "\tResolution perp. to the wire (mm): " << m_xy_resolution.value() << "\n";
  io << "\tCreate debug histograms: " << (m_create_debug_histos.value() ? "true" : "false") << "\n";
  if (true == m_create_debug_histos.value())
    io << "\t\t|--Name of output file with debug histograms: " << m_out_debug_filename.value() << "\n";

  return;
}


///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       CalculateNClusters       ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

bool DCHdigi_v01::IsParticleCreatedInsideDriftChamber(const edm4hep::MCParticle& thisParticle) const {
  auto vertex = thisParticle.getVertex(); // in mm
  auto vertexRsquared = vertex[0] * vertex[0] + vertex[1] * vertex[1];
  auto vertexZabs = std::fabs(vertex[2]);
  float DCH_halflengh = dch_data->Lhalf / dd4hep::mm;                // 2000;
  float DCH_rin_squared = std::pow(dch_data->rin / dd4hep::mm, 2);   // 350 * 350;
  float DCH_rout_squared = std::pow(dch_data->rout / dd4hep::mm, 2); // 2000 * 2000;
  return (vertexZabs < DCH_halflengh) && (vertexRsquared > DCH_rin_squared) && (vertexRsquared < DCH_rout_squared);
}

std::pair<uint32_t, std::vector<int> > DCHdigi_v01::CalculateClusters(const edm4hep::SimTrackerHit& input_sim_hit, TRandom3 & myRandom) const {

  const edm4hep::MCParticle& thisParticle = input_sim_hit.getParticle();
  // if gamma, optical photon, or other particle with null mass, or hit with zero energy deposited, return zero clusters
  if (22 == abs(thisParticle.getPDG()) || 0 == thisParticle.getMass() || 0 == input_sim_hit.getEDep())
    return {0., std::vector<int>{}};

  /// vector to accumulate the size (number of electrons) of each cluster
  std::vector<int> ClSz_vector;
  //_________________SET NECESSARY PARAMETERS FOR THE CLS ALGORITHM-----WALAA_________________//

  // Parameters needed for the CL Algo //Added by Walaa///////
  /// from Createclusters.cpp////
  float Eloss = 0.0;
  float EIzp = 15.8;
  float ExECl1 = 0;
  float cut = 1000; // controlla
  float EIzs = 25.6;
  float ExECl1totRec = 0;
  float rndCorr(0);
  const int nhEp = 10;
  float hEpcut[10] = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
  int minE = 1000;
  int maxE = 10000;
  int binE = 1000;
  int nhE = (maxE - minE) / binE;
  float maxEx0len(0), ExSgmlen(0);
  std::vector<double> CorrpMean;
  std::vector<double> CorrpSgm;
  std::vector<double> Corrdgmean;
  std::vector<double> Corrdgsgm;
  std::vector<double> Corrdglfrac;
  std::vector<double> Corrdglmpvl;
  std::vector<double> Corrdglsgml;
  std::vector<double> Corrdglmeang;
  std::vector<double> Corrdglsgmg;
  float maxEx0(0), maxExSlp(0), ExSgmlep(0), ExSgmhad(0);
  float MPVEx(0), SgmEx(0), MeanEx1(0), SgmEx1(0), frac(0), Slp(0), CorrSlp(0), CorrInt(0);

  /*________________________________________________________________________________*/
  bool IsSecondaryWithinDCH = IsParticleCreatedInsideDriftChamber(thisParticle);

  // Momentum from EDM4hep, in GeV
  double Momentum = sqrt((input_sim_hit.getMomentum().x * input_sim_hit.getMomentum().x) +
                         (input_sim_hit.getMomentum().y * input_sim_hit.getMomentum().y) +
                         (input_sim_hit.getMomentum().z * input_sim_hit.getMomentum().z));

  int thisparticle_pdgid = thisParticle.getPDG();
  int electron_pdgid = 11;
  constexpr double me = 0.511;         // electron mass in MeV
  constexpr double me_GeV = me * 1e-3; // electron mass in GeV

  /// number of clusters, with size > 1 electron
  int NClp(0);

  /// number of clusters, with size = 1 electron
  int NCl1(0);

  //___________________________________________________________________
  
  double thisparticle_mass = (thisParticle.getMass() / 1000.); // mass in GeV, required in MeV
  double bg = Momentum / thisparticle_mass;

  CorrpMean = flData->get_ClSzCorrpmean(bg);
  CorrpSgm = flData->get_ClSzCorrpsgm(bg);
  Corrdgmean = flData->get_ClSzCorrdgmean(bg);
  Corrdgsgm = flData->get_ClSzCorrdgsgm(bg);
  Corrdglfrac = flData->get_ClSzCorrdglfrac(bg);
  Corrdglmpvl = flData->get_ClSzCorrdglmpvl(bg);
  Corrdglsgml = flData->get_ClSzCorrdglsgml(bg);
  Corrdglmeang = flData->get_ClSzCorrdglmeang(bg);
  Corrdglsgmg = flData->get_ClSzCorrdglsgmg(bg);
  maxEx0 = flData->get_maxEx0(bg); // 379.4 electron 100 gev
  maxExSlp = flData->get_maxExSlp();

  MPVEx = flData->get_MPVExtra(bg);     // 103.2
  SgmEx = flData->get_SgmExtra(bg);     // 28.5;//
  MeanEx1 = flData->get_MeanExtra1(bg); // 13.8;//
  SgmEx1 = flData->get_SgmExtra1(bg);   // 9.72;//
  frac = flData->get_FracExtra1(bg);    // 0.13;//
  Slp = flData->get_SlopeExtra1(bg);    // 7.95;//
  CorrSlp = flData->get_ClSzCorrSlp(bg);
  CorrInt = flData->get_ClSzCorrInt(bg);

  double Tmax = (2.0 * me * pow(bg, 2) /
                 (1 + (2.0 * (1 + pow(bg, 2)) * me / thisparticle_mass) + pow(me / thisparticle_mass, 2))) *
                1e+6;
  float maxEcut = cut;
  if (Tmax < maxEcut) {
    maxEcut = Tmax;
  }
  //___________________________________________________________________
  if (not IsSecondaryWithinDCH) {
    // lepton PDG id goes from 11 to 16 (antiparticles have negative id)
    bool IsLepton = (11 <= abs(thisparticle_pdgid)) && (16 >= abs(thisparticle_pdgid));
    if (IsLepton) {
      ExSgmlep = flData->get_ExSgmlep();
    } else {
      // TODO Alvaro: gamma rays treated as hadrons, is that ok?
      ExSgmhad = flData->get_ExSgmhad();
    }

    /*________________________________________________________________________________*/

    TF1* land = new TF1("land", "landaun");
    land->SetParameter(0, 1);
    land->SetParameter(1, MPVEx);
    land->SetParameter(2, SgmEx);
    land->SetRange(0, maxEcut);

    TF1* exGauss = new TF1("exGauss", "[0]*([1]*TMath::Gaus(x,[2],[3],true)+(1.0-[1])*TMath::Exp(-x/[4])/[4])");
    exGauss->SetParameter(0, 1);
    exGauss->SetParameter(1, frac);
    exGauss->SetParameter(2, MeanEx1);
    exGauss->SetParameter(3, SgmEx1);
    exGauss->SetParameter(4, Slp);
    exGauss->SetRange(0, 90);

    float totExECl = 0.0;
    float ExECl = 0.0;
    double LengthTrack = input_sim_hit.getPathLength();
    double Etot_per_track = input_sim_hit.getEDep();
    LengthTrack *= 0.1;               // from mm to cm
    Eloss = Etot_per_track * 1.e09;   // from GeV to eV
    maxEx0len = maxEx0 * LengthTrack; // maxEx0 is a parameter
    if (IsLepton) {
      ExSgmlen = ExSgmlep * TMath::Sqrt(LengthTrack);
    } else {
      ExSgmlen = ExSgmhad * TMath::Sqrt(LengthTrack);
    }
    float maxExECl = (Eloss - maxEx0len + myRandom.Gaus(0, ExSgmlen)) / maxExSlp;

    if (maxExECl < EIzs) {
      maxExECl = 0.0;
    } // EIzs const = 25.6

    debug() << "Eloss= " << Eloss << "EIzs= " << EIzs << "EIzp= " << EIzp << "maxExECl= " << maxExECl << "totExECl"
            << totExECl << endmsg;

    // The following loop calculate number of clusters of size > 1, NClp
    for (int while1counter = 0; Eloss > (EIzp + EIzs) && maxExECl > totExECl && while1counter < 1e6; while1counter++) {
      ExECl = land->GetRandom(0, maxEcut, &myRandom);

      if (ExECl > EIzs) {
        Eloss -= EIzp;
        if (ExECl > (maxExECl - totExECl)) {
          ExECl = maxExECl - totExECl;
        }
        if (ExECl > Eloss)
          ExECl = Eloss;
        totExECl += ExECl;
        Eloss -= ExECl;
        NClp++;
        // CLSZ
        float tmpCorr = 0.0;
        for (int i = 0; i < nhEp; ++i) {
          if (ExECl >= (i == 0 ? 0 : hEpcut[i - 1]) && ExECl < hEpcut[i]) {
            tmpCorr = myRandom.Gaus(CorrpMean[i], CorrpSgm[i]);
          }
        }
        int ClSz = TMath::Nint(ExECl * CorrSlp + CorrInt - tmpCorr);
        if (ClSz < 2) {
          ClSz = 2;
        }
        ClSz_vector.push_back(ClSz);
      }
    }

    debug() << "Eloss= " << Eloss << "EIzp= " << EIzp << endmsg;

    // The following loop calculate number of clusters of size 1, NCl1
    for (int while2counter = 0; Eloss >= EIzp && while2counter < 1e6; while2counter++) {
      Eloss -= EIzp;
      ExECl1 = exGauss->GetRandom(&myRandom);
      if (ExECl1 > Eloss) {
        ExECl1 = Eloss;
      }
      ExECl1totRec += ExECl1;
      NCl1++;
      Eloss -= ExECl1;
      // cluster size is 1
      ClSz_vector.push_back(1);
    }

  } //-- end if particle is not secondary

  // if particle is a delta electron created inside the drift chamber
  int NCld(0);
  if (IsSecondaryWithinDCH && electron_pdgid == thisparticle_pdgid) {
    // 1 delta ray cause 1 cluster NCld (d=delta)
    NCld = 1;
    // Ekdelta in keV
    float Ekdelta = (TMath::Sqrt(Momentum * Momentum + me_GeV * me_GeV) - me_GeV) * 1e6;

    // TODO Alvaro: what is this for?
    {
      float tmpCl;
      int tmphE = (Ekdelta - minE) / binE;
      if (tmphE >= nhE)
        tmphE = nhE - 1;
      if (tmphE == nhE - 1) {
        rndCorr = myRandom.Uniform(0, 1);
        if (rndCorr < Corrdglfrac[tmphE]) {
          tmpCl = myRandom.Gaus(Corrdglmeang[tmphE], Corrdglsgmg[tmphE]);
        } else {
          tmpCl = myRandom.Landau(Corrdglmpvl[tmphE], Corrdglsgml[tmphE]);
        }
      } else {
        tmpCl = myRandom.Gaus(Corrdgmean[tmphE], Corrdgsgm[tmphE]);
      }

      int ClSz = TMath::Nint(Ekdelta * CorrSlp + CorrInt - tmpCl);
      // TODO Alvaro: should it be 1 instead?
      if (ClSz < 2) {
        ClSz = 2;
      }
      ClSz_vector.push_back(ClSz);
    }
  }
  // conclusion: if hit caused by delta electron, NCld = 1, number of electrons ClSz=2

  // value to be returned, sum of number of clusters size=1, size>1, and clusters caused by delta rays
  int total_number_of_clusters = NCl1 + NClp + NCld;
  debug() << "Ncl= " << total_number_of_clusters << " NCl1= " << NCl1 << "NClp= " << NClp << "NCld= " << NCld << endmsg;

  if (ClSz_vector.size() != std::size_t(total_number_of_clusters))
    debug() << "Array of cluster sizes does not match total number of clusters\n";

  // return {total_number_of_clusters, total_number_of_electrons_over_all_clusters};
/*  
  //---------changes made by Muhammad saiel---------//
  // ---- [DEBUG: Walaa-based cluster stats] ----
if (m_create_debug_histos)
{
    const double l_mm = input_sim_hit.getPathLength();
    const double l_cm = 0.1 * l_mm;
    const double meanSpacing_mm = (total_number_of_clusters > 0)
                                  ? (l_mm / total_number_of_clusters)
                                  : 0.0;

    hNcl_perStep_Walaa->Fill(total_number_of_clusters);
    hNcl_vs_l_Walaa->Fill(l_cm, total_number_of_clusters);
    if (meanSpacing_mm > 0.0) hClSpacing_Walaa_mm->Fill(meanSpacing_mm);
}
//------------------End of changes---------//
*/  
  return {total_number_of_clusters, ClSz_vector};
}
