
#include "WireTrackerDigi_v02.h"
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

// Gaudi and DD4hep
#include <GaudiKernel/MsgStream.h>
#include "DD4hep/Detector.h"
#include "DD4hep/DetElement.h"

///////////////////////////////////////////////////////////////////////////////////////
////////////////////////////        constructor       /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
WireTrackerDigi_v02::WireTrackerDigi_v02(const std::string& name, ISvcLocator* svcLoc)
  :MultiTransformer(name, svcLoc,
    {
      //input collection
      KeyValues("InputClusterCollection", {"GasIonizationClusters"}),
      KeyValues("HeaderName", {"EventHeader"}),
    },
    {
      //output collection
		  KeyValues("OutputDigiHitCollection", {"SenseWire_DigiCollection"}),
			KeyValues("DigiSimAssociationCollection", {"SenseWire_DigiSimAssociationCollection"}),
    } 
  )  
  {
    m_geoSvc = serviceLocator()->service(m_geoSvcName);	// Access detector geometry service
    m_uidSvc = serviceLocator()->service(m_uidSvcName);	// Access reproducible unique-ID seed service  
  }
//--------------------------------------------------------------------------------------------------//

// Build an event-reproducible random-number generator from the EventHeader.
TRandom3 WireTrackerDigi_v02::CreateRandomEngine(
		const edm4hep::EventHeaderCollection& headers) const {
  const auto seed = m_uidSvc->getUniqueID(headers, this->name());
  TRandom3 rng(seed);
  return rng;
}
//------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       initialize       ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
StatusCode WireTrackerDigi_v02::initialize() 
{
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

  // Retrieve the detector element corresponding to the subdetector name.
  dd4hep::DetElement detectorElement = m_geoSvc->getDetector()->detectors().at(detectorName);

	// Retrieve the wire-tracker extension using the DD4hep base extension key.
	auto* wireTracker_info = detectorElement.extension<dd4hep::rec::WireTracker_info_struct>();
	if (!wireTracker_info) {
		error() << "No WireTracker_info_struct extension was attached to detector <<"
			<< detectorName << ">>." << endmsg;
		return StatusCode::FAILURE;
	}

	// Cast the generic wire-tracker extension to the concrete 
  // geometry helper used by this detector.
	dch_data = dynamic_cast<dd4hep::rec::DCH_info*>(wireTracker_info);
	if (!dch_data || !dch_data->IsValid()) {
		error() << "No valid wire-tracker geometry extension was found for detector <<"
				<< detectorName << ">>." << endmsg;
		return StatusCode::FAILURE;
	}

  // Retrieve the readout associated with the detector element (subdetector)
  dd4hep::SensitiveDetector sensitiveDetector = m_geoSvc->getDetector()->sensitiveDetector(detectorName);
  if (!sensitiveDetector.isValid()) {
    error() << "No valid Sensitive Detector was found for detector <<" << detectorName << ">>." << endmsg;
    return StatusCode::FAILURE;
  }

  // Retrieve the readout associated with the sensitive detector.
  dd4hep::Readout readout = sensitiveDetector.readout();

  // Retrieve the CellID decoder from the readout.
  m_decoder = readout.idSpec().decoder();
  if (!m_decoder) {
  error() << "Failed to retrieve CellID decoder from detector readout <<" << readout.name() << ">>." << endmsg;
  return StatusCode::FAILURE;
  }

  // Check that the resolution parameters are non-negative.
  if (0 > m_z_resolution_mm.value()) {
    error() << "Z resolution input value can not be negative!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (0 > m_xy_resolution_mm.value()) {
    error() << "Radial (XY) resolution input value can not be negative!" << endmsg;
    return StatusCode::FAILURE;
  }

  ///=================================================================================//
  //                            initialize x-t relation                               //
  //==================================================================================// 
  // Load the 2D x-t lookup table used to convert cluster position into drift time.
  try {
	  m_xt2d = std::make_unique<DCHXT2DLUT>(m_xtFileName.value(), "xt_mean", "xt_error");
  } catch (const std::exception& e) {
	  error() << e.what() << endmsg;
	  return StatusCode::FAILURE;
  }
  info() << "[XT2D] Loaded XT LUT from " << m_xtFileName.value() << endmsg;

  // Print the configuration to the log.
  std::stringstream ss;
  PrintConfiguration(ss);
  info() << ss.str().c_str() << endmsg;

  //===================================================================================//
  //                           Create histograms for debugging:                        // 
  //===================================================================================//
  if (m_create_debug_histos.value()) 
  {
    m_histSvc = service("THistSvc");
    if (!m_histSvc) {
      error() << "Unable to locate Histogram Service" << endmsg;
      return StatusCode::FAILURE;
    }

    hNe = new TH1F("hNe", "Cluster Size Distribution;Cluster size;Entries",15, 0.0, 15.0);
    if (m_histSvc->regHist("/rec/hNe_clusters", hNe).isFailure()) {
      error() << "Couldn't register hNe histogram" << endmsg;
      return StatusCode::FAILURE;
    }
    hTd = new TH1F  ("hTd", "Drift Time; t (ns); Entries", 50, 0, 1000);
    if (m_histSvc->regHist("/rec/hTd", hTd).isFailure()) {
      error() << "Couldn't register hTd histogram" << endmsg;
      return StatusCode::FAILURE;
    }
    hXT = new TProfile ("hXT", "x-t relation; r (cm); t (ns)", 50, 0, 1.0);
    if (m_histSvc->regHist("/rec/hXT", hXT).isFailure()) {
      error() << "Couldn't register hXT histogram" << endmsg;
      return StatusCode::FAILURE;
    }
    hXYT = new TProfile2D("hXYT","Drift time;X (cm);Y (cm);Time (ns)",
      50, -0.9, 0.9, 50, -0.9, 0.9);
    if (m_histSvc->regHist("/rec/hXYT", hXYT).isFailure()) {
      error() << "Couldn't register hXYT histogram" << endmsg;
      return StatusCode::FAILURE;
    }

  } //End of histogram definition

  return StatusCode::SUCCESS;
} //end of Initialize

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       operator()       ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
std::tuple<
	edm4hep::SenseWireHitCollection,
	edm4hep::TrackerHitSimTrackerHitLinkCollection>
WireTrackerDigi_v02::operator()(const edm4hep::SimTrackerHitCollection& input_cluster_hits,
		const edm4hep::EventHeaderCollection&   headers) const 
{
  // Reproducible event-level random engine
  auto myRandom = this->CreateRandomEngine(headers);

  // Create the collections we are going to return
  edm4hep::SenseWireHitCollection output_digi_hits;
  edm4hep::TrackerHitSimTrackerHitLinkCollection output_digi_sim_association;

  // Per-cell information needed to create one final SenseWireHit.
  struct CellAccum {
    std::vector<double> times_ns;	// Absolute electron arrival times
    std::vector<uint16_t> electronsPerCluster;		// number of electron per cluster
    std::vector<edm4hep::SimTrackerHit> simHits;	// clusters linked to the final cell hit

    // Representative cluster metadata for the final cell-level hit.
    // Policy: use the earliest cluster time seen in this cell.
    bool hasRepresentative = false;
    double representativeClusterTime_ns = std::numeric_limits<double>::max();

    // cluster position and distance from the nearest sense wire
    edm4hep::Vector3d positionSW{};
    float distanceToWire = 0.f;
    float wireAzimuthalAngle = 0.f;
    float wireStereoAngle = 0.f; 

    float eDep = 0.f;
  };
  std::unordered_map<uint64_t, CellAccum> cellAccumMap;

  // loop over hit collection
  // Each input object represents one cluster
  int loop_index = 0;
  for (const auto& input_cluster_hit : input_cluster_hits) 
  {
    loop_index++;

    // Get the cellID of the cluster hit.
    dd4hep::DDSegmentation::CellID cellid = input_cluster_hit.getCellID();

    // Decode the cellID to get layer and cell indices.
    const int superlayer = this->CalculateSuperLayerFromCellID(cellid);
    const int sector = this->CalculateSectorFromCellID(cellid);
    const int ilayer = this->CalculateLayerFromCellID(cellid);
    const int nphi = this->CalculateNphiFromCellID(cellid);

    // Get the cluster position in global coordinates (cm).
    //auto hit_position = Convert_EDM4hepVector_to_TVector3(input_cluster_hit.getPosition(),MM_TO_CM);

    auto hit_position = Convert_EDM4hepVector_to_XYZVector(input_cluster_hit.getPosition(), MM_TO_CM);
    auto hit_to_wire_vector = dch_data->Calculate_hitpos_to_wire_vector(
      superlayer, ilayer, sector, nphi, hit_position
    );
    auto hit_projection_on_the_wire = hit_position + hit_to_wire_vector;
    auto wire_direction_ez = dch_data->Calculate_wire_vector_ez(superlayer, ilayer, sector, nphi);

    // Debug printout for the clusters
    debug() 
      << "Cluster no." <<loop_index 
      <<" Layer = " << ilayer 
      <<" cell = " << nphi 
      << endmsg;

    // Sample the number of electrons associated with this cluster.
    const int cluster_size = sample_cluster_size(myRandom);
    if (cluster_size <= 0) {
      continue;
    }
    if (m_create_debug_histos.value() && hNe) {
	    // Fill debug histogram
	    hNe->Fill(cluster_size);
    }
    
    //===================================================================================//
    //                         Hit to wire distance Calculation                          //
    //===================================================================================//
    // Project the cluster position onto the corresponding sense wire.
    // smear position along the wire
    double smearing_z = myRandom.Gaus(0.0, m_z_resolution_mm.value() * MM_TO_CM);

    // Add smearing along the wire direction
    hit_projection_on_the_wire += smearing_z * (wire_direction_ez.Unit());

    // Smear the transverse distance to the sense wire.
    double smearing_xy = myRandom.Gaus(0.0, m_xy_resolution_mm.value() * MM_TO_CM);
    const double distanceToWire_Magnitude = std::sqrt(hit_to_wire_vector.Dot(hit_to_wire_vector));

    // protect against negative values
    const double distanceToWire_smeared = std::max(0.0, distanceToWire_Magnitude + smearing_xy);
    
    // Convert the smeared wire position back to EDM4hep units [mm].
    auto positionSW = Convert_XYZVector_to_EDM4hepVector(hit_projection_on_the_wire, 1.0 / MM_TO_CM);
    const float distanceToWire = distanceToWire_smeared / MM_TO_CM;

    // Store the wire direction parameters required by edm4hep::SenseWireHit.
    // Together, wireAzimuthalAngle and wireStereoAngle define the local wire axis.
    const auto WireAzimuthalAngle = this->dch_data->Get_cell_phi_angle(superlayer, ilayer, sector, nphi);
    float WireStereoAngle = 0;
    {
      auto l = this->dch_data->database.at(ilayer);
      // Radius of the cell center at z = 0
      auto cell_rave_z0 = 0.5 * (l.radius_fdw_z0 + l.radius_fuw_z0);
      // DD4hep twisted-tube convention gives the wire stereo angle with opposite sign.
      WireStereoAngle = (-1.) * dch_data->StereoSign(l) * dch_data->stereoangle_z0(cell_rave_z0);
    }

    //===================================================================================//
    //                                 Drift time  Calculation                           //
    //===================================================================================//
    // Store per-electron drift time
    std::vector<double> electron_times_ns;
    electron_times_ns.reserve(cluster_size);

    auto cluster_pos = hit_position;

    auto cl_to_wire = dch_data->Calculate_hitpos_to_wire_vector(
        superlayer, ilayer, sector, nphi, cluster_pos);

    auto u = wire_direction_ez.Unit();

    // Build an orthonormal basis in the plane transverse to the sense wire.
    ROOT::Math::XYZVector e1(1.0, 0.0, 0.0);
    e1 -= e1.Dot(u) * u;

    if (std::sqrt(e1.Dot(e1)) < 1e-9) {
      e1 = ROOT::Math::XYZVector(0.0, 1.0, 0.0);
      e1 -= e1.Dot(u) * u;
    }

    e1 = e1.Unit();

    auto e2 = u.Cross(e1).Unit();

    auto v_perp = cl_to_wire - cl_to_wire.Dot(u) * u;

    // Local transverse coordinates used to query the x-t lookup table.
    const double x_cluster_cm = v_perp.Dot(e1);
    const double y_cluster_cm = v_perp.Dot(e2);
    const double r_cluster_cm = std::sqrt(v_perp.Dot(v_perp));

    // Sample one drift time per electron in the cluster
    for (int ie = 0; ie < cluster_size; ++ie) 
    {
      const double t_ns = m_xt2d->sampleTimeNs(x_cluster_cm, y_cluster_cm, myRandom);
      electron_times_ns.push_back(t_ns);
      if (m_create_debug_histos.value()) {
        hTd->Fill(t_ns);
        hXT->Fill(r_cluster_cm, t_ns);
        hXYT->Fill(x_cluster_cm, y_cluster_cm, t_ns);
      }
    }

    //=====================================================================================//
    // 			                  Accumulate cluster information per cellID				             //
    //=====================================================================================//			
    auto& cellAccum = cellAccumMap[static_cast<uint64_t>(cellid)];

    // Keep the cluster-level SimTrackerHit contributing to this final cell hit.
    const auto alreadyStored = std::find_if(
      cellAccum.simHits.begin(),
      cellAccum.simHits.end(),
      [&](const edm4hep::SimTrackerHit& h) {
        return h.getObjectID().index == input_cluster_hit.getObjectID().index &&
        h.getObjectID().collectionID == input_cluster_hit.getObjectID().collectionID;
      }
    );
    if (alreadyStored == cellAccum.simHits.end()) {
      cellAccum.simHits.push_back(input_cluster_hit);
    }
  
    // Absolute time of this cluster
    const double clusterTime_ns = input_cluster_hit.getTime();
    
    // Keep total deposited energy summed over all clusters contributing to this cell.
    cellAccum.eDep += static_cast<float>(input_cluster_hit.getEDep());

    // Store representative DigiHit metadata once per cell.
    if (!cellAccum.hasRepresentative || clusterTime_ns < cellAccum.representativeClusterTime_ns) {

      // Use the earliest cluster as representative for position and geometry metadata.
      cellAccum.hasRepresentative = true;
      cellAccum.representativeClusterTime_ns = clusterTime_ns;
      cellAccum.positionSW = positionSW;
      cellAccum.distanceToWire = distanceToWire;
      cellAccum.wireAzimuthalAngle = WireAzimuthalAngle;
      cellAccum.wireStereoAngle = WireStereoAngle;
    } 
    // Store absolute electron times in the event time frame
    for (double tdrift_ns : electron_times_ns) {
      cellAccum.times_ns.push_back(clusterTime_ns + tdrift_ns);
    }
    // One input gas cluster contributes one cluster-size entry.  
    cellAccum.electronsPerCluster.push_back(static_cast<uint16_t>(cluster_size));
    
  } // end loop over hit collection

  //====================================================================================//
  //                    Create one final SenseWireHit per cell                          //        
  //====================================================================================//
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

    // Final DigiHit time is the earliest drifted-electron arrival in this cell.
    const auto it_min = std::min_element(acc.times_ns.begin(), acc.times_ns.end());  
    const double t_hit_ns = *it_min;  
    auto digiHit = output_digi_hits.create();
      
    digiHit.setCellID(kv.first);
    digiHit.setType(0);
    digiHit.setQuality(0);
  
    digiHit.setTime(t_hit_ns);
    digiHit.setEDep(acc.eDep);
    digiHit.setEDepError(0.0f);
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
  return std::make_tuple( std::move(output_digi_hits), 
    std::move(output_digi_sim_association));
  
} // End of the operator Block

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       finalize       //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
StatusCode WireTrackerDigi_v02::finalize() 
{
  return StatusCode::SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////       PrintConfiguration       ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
void WireTrackerDigi_v02::PrintConfiguration(std::ostream& io) const
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
