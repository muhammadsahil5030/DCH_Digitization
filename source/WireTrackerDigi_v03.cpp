
#include "WireTrackerDigi_v03.h"
#include "DCHUtilities.h"
#include "DCHXT2DLUT.h"

#include "TMath.h"
#include "TProfile.h"

// STL
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <sstream>
#include <unordered_map>

// Gaudi
#include <GaudiKernel/MsgStream.h>

#include "DD4hep/Detector.h"
#include "DD4hep/DetElement.h"


///////////////////////////////////////////////////////////////////////////////////////
////////////////////////////        constructor       /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
WireTrackerDigi_v03::WireTrackerDigi_v03(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc,
       {
        //input collection
        KeyValues("InputClusterCollection", {"GasIonizationClusters"}),
        KeyValues("HeaderName", {"EventHeader"}),
        },{
          //output collection
		      KeyValues("OutputDigiHitCollection", {"SenseWire_DigiCollection"}),
			    KeyValues("DigiSimAssociationCollection", {"SenseWire_DigiSimAssociationCollection"}),
        })
        {
          m_geoSvc = serviceLocator()->service(m_geoSvcName);	// Access detector geometry service
          m_uidSvc = serviceLocator()->service(m_uidSvcName);	// Access reproducible unique-ID seed service  
        }

// Random number generator
TRandom3 WireTrackerDigi_v03::CreateRandomEngine(
		const edm4hep::EventHeaderCollection& headers) const {
  const auto seed = m_uidSvc->getUniqueID(headers, this->name());
  TRandom3 rng(seed);
  return rng;
}
//------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       initialize       ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
StatusCode WireTrackerDigi_v03::initialize() 
{
  if (!m_uidSvc) {
    error() << "Unable to get UniqueIDGenSvc" << endmsg;
    return StatusCode::FAILURE;
  }
  if (!m_geoSvc) {
    error() << "Unable to get GeoSvc" << endmsg;
    return StatusCode::FAILURE;
  }

  if (0 > m_z_resolution_mm.value()) {
    error() << "Z resolution input value can not be negative!" << endmsg;
    return StatusCode::FAILURE;
  }

  if (0 > m_xy_resolution_mm.value()) {
    error() << "Radial (XY) resolution input value can not be negative!" << endmsg;
    return StatusCode::FAILURE;
  }

  // Retrieve the subdetector
  std::string detectorName(m_detectorName.value());
  if (0 == m_geoSvc->getDetector()->detectors().count(detectorName)) {
    error() << "Detector <<" << detectorName <<  ">> does not exist." << endmsg;
    return StatusCode::FAILURE;
  }

  dd4hep::DetElement detectorElement = m_geoSvc->getDetector()->detectors().at(detectorName);

  // Retrieve the DCH geometry helper extension
  this->dch_data = detectorElement.extension<dd4hep::rec::DCH_info>();
  if (!dch_data || !dch_data->IsValid()) {
    error() << "No valid data extension was found for detector <<" << detectorName << ">>." << endmsg;
    return StatusCode::FAILURE;
  }
	
  // Retrieve the readout associated with the detector element (subdetector)
  dd4hep::SensitiveDetector sensitiveDetector = m_geoSvc->getDetector()->sensitiveDetector(detectorName);
  if (!sensitiveDetector.isValid()) {
    error() << "No valid Sensitive Detector was found for detector <<" << detectorName << ">>." << endmsg;
    return StatusCode::FAILURE;
  }

  dd4hep::Readout readout = sensitiveDetector.readout();
  // set the cellID decoder
  m_decoder = readout.idSpec().decoder();
  if (!m_decoder) {
  error() << "Failed to retrieve CellID decoder from detector readout <<" << readout.name() << ">>." << endmsg;
  return StatusCode::FAILURE;
  }

  ///////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////  initialize x-t relation  //////////////////////////
  /////////////////////////////////////////////////////////////////////////////////// 
  // Load the 2D x-t lookup table used to convert cluster position into drift time.

  try {
	  m_xt2d = std::make_unique<DCHXT2DLUT>(m_xtFileName.value(), "xt_mean", "xt_error");
  } catch (const std::exception& e) {
	  error() << e.what() << endmsg;
	  return StatusCode::FAILURE;
  }
  info() << "[XT2D] Loaded XT LUT from " << m_xtFileName.value() << endmsg;
  //-----------------------------------------------------------------------------------//

  std::stringstream ss;
  PrintConfiguration(ss);
  info() << ss.str().c_str() << endmsg;

  // Create histograms for debugging:
  if (m_create_debug_histos.value()) 
  {
    m_histSvc = service("THistSvc");
    if (!m_histSvc) {
      error() << "Unable to locate Histogram Service" << endmsg;
      return StatusCode::FAILURE;
    }

    // Cluster size distribution (number of electrons per cluster)
    hNe = new TH1F("hNe", "Cluster Size Distribution;Cluster size;Entries",15, 0.0, 15.0);
    if (m_histSvc->regHist("/rec/hNe_clusters", hNe).isFailure()) {
      error() << "Couldn't register hNe histogram" << endmsg;
      return StatusCode::FAILURE;
    }

    // Drift time 
    hTd = new TH1F  ("hTd", "Drift Time; t (ns); Entries", 50, 0, 1000);
    if (m_histSvc->regHist("/rec/hTd", hTd).isFailure()) {
      error() << "Couldn't register hTd histogram" << endmsg;
      return StatusCode::FAILURE;
    }

    // Drift time vs radius
    hXT = new TProfile ("hXT", "x-t relation; r (cm); t (ns)", 50, 0, 1.0);
    if (m_histSvc->regHist("/rec/hXT", hXT).isFailure()) {
      error() << "Couldn't register hXT histogram" << endmsg;
      return StatusCode::FAILURE;
    }

    // Drift time in XY plane
    hXYT = new TProfile2D("hXYT","Drift time;X (cm);Y (cm);Time (ns)",
      50, -0.9, 0.9, 50, -0.9, 0.9);
    if (m_histSvc->regHist("/rec/hXYT", hXYT).isFailure()) {
      error() << "Couldn't register hXYT histogram" << endmsg;
      return StatusCode::FAILURE;
    }

} //End of histogram definition
//----------------------------------------------------------------------------------------//

  return StatusCode::SUCCESS;

} //end of Initialize

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       operator()       ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

std::tuple<
	edm4hep::SenseWireHitCollection,
	edm4hep::TrackerHitSimTrackerHitLinkCollection>
WireTrackerDigi_v03::operator()(const edm4hep::SimTrackerHitCollection& input_cluster_hits,
		const edm4hep::EventHeaderCollection&   headers) const 
{
  /// initialize engines
  auto myRandom = this->CreateRandomEngine(headers);

  // Create the collections we are going to return
  edm4hep::SenseWireHitCollection output_digi_hits;
  edm4hep::TrackerHitSimTrackerHitLinkCollection output_digi_sim_association;

  //Varible initailizations:
  int ilayer = -1;
  int nphi = -1;

  // loop over hit collection
  int loop_index = 0;

  // Per-cell information needed to create one final SenseWireHit.
  struct CellAccum {
    std::vector<double> times_ns;	// Absolute electron arrival times
    std::vector<uint16_t> electronsPerCluster;		// number of electron per cluster
    std::vector<edm4hep::SimTrackerHit> simHits;	// Keep all contributing sim hits

    // Representative cluster metadata for the final cell-level hit.
    // Policy: use the earliest cluster time seen in this cell.
    bool hasRepresentative = false;
    double representativeClusterTime_ns = std::numeric_limits<double>::max();

    edm4hep::Vector3d positionSW{};
    float distanceToWire = 0.f;
    float wireAzimuthalAngle = 0.f;
    float wireStereoAngle = 0.f;

    std::int32_t type = 0;
    std::int32_t quality = 0;

    float eDep = 0.f;
    float eDepError = 0.f;
  };
  std::unordered_map<uint64_t, CellAccum> cellAccumMap;

  // Each input object represents one generated gas-ionization cluster
  for (const auto& input_cluster_hit : input_cluster_hits) 
  {
    loop_index++;

    dd4hep::DDSegmentation::CellID cellid = input_cluster_hit.getCellID();

    ilayer = this->CalculateLayerFromCellID(cellid);
    nphi = this->CalculateNphiFromCellID(cellid);

    auto hit_position = Convert_EDM4hepVector_to_TVector3(input_cluster_hit.getPosition(),MM_TO_CM);

        debug() << "Cluster no." <<loop_index
                <<" Layer = " << ilayer
                <<" cell = " << nphi
                << endmsg;

    // Sample the number of electrons associated with this gas cluster.
    const int Ne_cluster = sample_cluster_size(myRandom);
    if (Ne_cluster <= 0) {
            continue;
    }
    if (m_create_debug_histos.value() && hNe) {
    
	    // Fill debug histogram
	    hNe->Fill(Ne_cluster);
    }
    //-------------------------------------------------------------------------
    
    // Project the cluster position onto the corresponding sense wire.
    TVector3 hit_to_wire_vector = this->dch_data->Calculate_hitpos_to_wire_vector(
		    ilayer, nphi, hit_position);

    TVector3 hit_projection_on_the_wire = hit_position + hit_to_wire_vector;
    if (m_create_debug_histos.value()) {
      double distance_hit_wire = hit_to_wire_vector.Mag();
      hDpw->Fill(distance_hit_wire);
    }
    TVector3 wire_direction_ez = this->dch_data->Calculate_wire_vector_ez(ilayer, nphi);
    // -------------------------------------------------------------------------
    // smear position along the wire
    double smearing_z = myRandom.Gaus(0.0, m_z_resolution_mm.value() * MM_TO_CM);
    if (m_create_debug_histos.value())
      hSz->Fill(smearing_z);

    hit_projection_on_the_wire += smearing_z * (wire_direction_ez.Unit());
    if (m_create_debug_histos.value()) {
      // the distance from the hit projection and the wire should be zero
      TVector3 dummy_vector = this->dch_data->Calculate_hitpos_to_wire_vector(
		      ilayer, nphi, hit_projection_on_the_wire);
      hDww->Fill(dummy_vector.Mag());
    }

    // smear position perpendicular to the wire
    double smearing_xy = myRandom.Gaus(0.0, m_xy_resolution_mm.value() * MM_TO_CM);
    if (m_create_debug_histos.value())
      hSxy->Fill(smearing_xy);
    const double distanceToWire_real = hit_to_wire_vector.Mag();

    // protect against negative values
    const double distanceToWire_smeared = std::max(0.0, distanceToWire_real + smearing_xy);

    const std::int32_t type = 0;
    const std::int32_t quality = 0;
    const float eDepError = 0;
    // length units back to mm
    auto positionSW = Convert_TVector3_to_EDM4hepVector(hit_projection_on_the_wire, 1. / MM_TO_CM);
    // auto  directionSW    = Convert_TVector3_to_EDM4hepVector(wire_direction_ez, 1. / MM_TO_CM);
    const float distanceToWire = distanceToWire_smeared / MM_TO_CM;

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
//                                      Drift time                                  //
//====================================================================================
    // Store per-electron drift time
    std::vector<double> electron_times_ns;
    electron_times_ns.reserve(Ne_cluster);

    // Cluster position is already provided by the upstream clusterizer
    TVector3 cluster_pos = hit_position;

    // Compute vector from cluster position to the closest point on the wire
    TVector3 cl_to_wire = this->dch_data->Calculate_hitpos_to_wire_vector(    
		    ilayer, nphi, cluster_pos);

    // Local 2D transverse basis for the x-t lookup
    TVector3 u = wire_direction_ez.Unit();
    TVector3 e1(1.0, 0.0, 0.0);
    e1 -= (e1.Dot(u)) * u;

    if (e1.Mag() < 1e-9) {
	    e1 = TVector3(0.0, 1.0, 0.0);
	    e1 -= (e1.Dot(u)) * u;
    }

    e1 = e1.Unit();
    TVector3 e2 = u.Cross(e1).Unit();

    // Perpendicular displacement from cluster to wire
    TVector3 v_perp = cl_to_wire - (cl_to_wire.Dot(u)) * u;

    // True 2D coordinates in the wire-centered transverse plane [cm]
    const double x_cluster_cm = v_perp.Dot(e1);
    const double y_cluster_cm = v_perp.Dot(e2);

    // Radius kept for monitoring
    const double r_cluster_cm = v_perp.Mag();

    // Sample one drift time per electron in the cluster
    for (int ie = 0; ie < Ne_cluster; ++ie) 
	{
    const double t_ns = m_xt2d->sampleTimeNs(x_cluster_cm, y_cluster_cm, myRandom);

    electron_times_ns.push_back(t_ns);

    if (m_create_debug_histos.value()) {
        hTd->Fill(t_ns);
        hXT->Fill(r_cluster_cm, t_ns);
        if (hXYT) hXYT->Fill(x_cluster_cm, y_cluster_cm, t_ns);
    }
}

//=======================================================================================
// 			                  Accumulate cluster information per cellID				             //
//======================================================================================= 			

    auto& cellAccum = cellAccumMap[static_cast<uint64_t>(cellid)];

    // Keep the cluster-level SimTrackerHit contributing to this final cell hit.
    const auto alreadyStored = std::find_if(
		    cellAccum.simHits.begin(),
		    cellAccum.simHits.end(),
		    [&](const edm4hep::SimTrackerHit& h) {
		    return h.getObjectID().index == input_cluster_hit.getObjectID().index &&
		    h.getObjectID().collectionID == input_cluster_hit.getObjectID().collectionID;
		    });
    if (alreadyStored == cellAccum.simHits.end()) {
	    cellAccum.simHits.push_back(input_cluster_hit);
    }

    // Absolute time of this cluster
    const double clusterTime_ns = input_cluster_hit.getTime();

    // Keep total deposited energy summed over all clusters contributing to this cell.
    cellAccum.eDep += static_cast<float>(input_cluster_hit.getEDep());

    // Store representative DigiHit metadata once per cell.
    if (!cellAccum.hasRepresentative ||
		    clusterTime_ns < cellAccum.representativeClusterTime_ns) {
	    cellAccum.hasRepresentative = true;
	    cellAccum.representativeClusterTime_ns = clusterTime_ns;
	    cellAccum.positionSW = positionSW;
	    cellAccum.distanceToWire = distanceToWire;
	    cellAccum.wireAzimuthalAngle = WireAzimuthalAngle;
	    cellAccum.wireStereoAngle = WireStereoAngle;
	    cellAccum.type = type;
	    cellAccum.quality = quality;
	    cellAccum.eDepError = eDepError;
    }

    // Store absolute electron times in the event time frame
    for (double tdrift_ns : electron_times_ns) {
	    cellAccum.times_ns.push_back(clusterTime_ns + tdrift_ns);
    }

    // One input gas cluster contributes one cluster-size entry.
    cellAccum.electronsPerCluster.push_back(static_cast<uint16_t>(Ne_cluster));
  
  } // end loop over hit collection

//====================================================================================
//                    Create one final SenseWireHit per cell                        //        
//====================================================================================

  // Retrieve accumulated data per cell (all electrons)
  for (const auto& kv : cellAccumMap) 
  {
  const auto& acc = kv.second;

  if (acc.times_ns.empty()) {
    continue;
  }

  if (!acc.hasRepresentative) {
    continue;
  }

  const auto it_min =
      std::min_element(acc.times_ns.begin(), acc.times_ns.end());

  const double t_hit_ns = *it_min;

    auto digiHit = output_digi_hits.create();

  digiHit.setCellID(kv.first);
  digiHit.setType(acc.type);
  digiHit.setQuality(acc.quality);

  digiHit.setTime(t_hit_ns);
  digiHit.setEDep(acc.eDep);
  digiHit.setEDepError(acc.eDepError);
  digiHit.setPosition(acc.positionSW);

  digiHit.setPositionAlongWireError(m_z_resolution_mm.value());
  digiHit.setWireAzimuthalAngle(acc.wireAzimuthalAngle);
  digiHit.setWireStereoAngle(acc.wireStereoAngle);
  digiHit.setDistanceToWire(acc.distanceToWire);
  digiHit.setDistanceToWireError(m_xy_resolution_mm.value());

  for (auto ne : acc.electronsPerCluster) {
    digiHit.addToNElectrons(ne);
  }

  for (const auto& simHit : acc.simHits) {
    auto link = output_digi_sim_association.create();
    link.setWeight(1.0f);
    link.setFrom(digiHit);
    link.setTo(simHit);
  }
} // End of loop over accumulated cells

// Return direct DigiHits and their truth links for this event.
  return std::make_tuple(
      std::move(output_digi_hits), 
      std::move(output_digi_sim_association));

} // End of the operator Block

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       finalize       //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
StatusCode WireTrackerDigi_v03::finalize() 
{
  return StatusCode::SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       PrintConfiguration       ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
void WireTrackerDigi_v03::PrintConfiguration(std::ostream& io) const
{
  io << "Direct wire-tracker DigiHit producer will use the following components:\n";
  io << "\tGeometry Service: " << m_geoSvcName.value().c_str() << "\n";
  io << "\tUID Service: " << m_uidSvcName.value().c_str() << "\n";
  io << "\tDetector name: " << m_detectorName.value().c_str() << "\n";
  io << "\t\t|--Volume bitfield: " << m_decoder->fieldDescription().c_str() << "\n";
  io << "\t\t|--Number of layers: " << dch_data->database.size() << "\n";
  io << "\tResolution along the wire (mm): " << m_z_resolution_mm.value() << "\n";
  io << "\tResolution perp. to the wire (mm): " << m_xy_resolution_mm.value() << "\n";
    io << "\tCreate debug histograms (THistSvc): " 
      << (m_create_debug_histos.value() ? "true" : "false") << "\n";
  return;
}
