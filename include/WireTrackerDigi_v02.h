/** ======= WireTrackerDigi_v02 ==========
 * Gaudi Algorithm for direct wire-tracker DigiHit creation
 *
 * @author Muhammad Saiel, Giovanni Francesco Tassielli, Nicola De Filippis
 * @date   2026-04
 *
 * <h4>Purpose</h4>
 * This algorithm is the direct DigiHit-producing stage of the split IDEA DCH /
 * wire-tracker digitization chain.
 *
 * It consumes gas-ionization clusters produced upstream by
 * DistributeClustersInGas and creates ready-to-use edm4hep::SenseWireHit
 * objects. Each output SenseWireHit represents the direct digitized hit
 * information for one readout cell.
 *
 * This algorithm does not create analog waveforms, electronics response,
 * electronic noise, or ADC digitized waveforms. Those detector-response
 * waveform steps are handled separately by WaveformFromWireTracker.
 *
 * This code uses DD4hep natural length units (cm), while EDM4hep positions and
 * distances are stored in mm. Please handle unit conversions carefully. <br>
 *
 * <br><h4>Input collections and prerequisites</h4>
 * Processor requires:
 *   - edm4hep::SimTrackerHitCollection produced upstream by
 *     DistributeClustersInGas, where each SimTrackerHit represents one
 *     gas-ionization cluster
 *   - edm4hep::EventHeaderCollection for reproducible random seeds
 *
 * <br><h4>Scope of this stage</h4>
 * This algorithm is responsible for:
 *   - reading gas-ionization clusters
 *   - projecting each cluster position onto the corresponding sense wire
 *   - applying position smearing
 *   - computing drift times from the x-t lookup table
 *   - sampling the electron multiplicity for each gas cluster
 *   - accumulating clusters/electrons belonging to the same cell
 *   - producing one final edm4hep::SenseWireHit per readout cell
 *   - producing edm4hep::TrackerHitSimTrackerHitLink truth links
 *
 * <br><h4>Output</h4>
 * Processor produces:
 *   - edm4hep::SenseWireHitCollection
 *   - edm4hep::TrackerHitSimTrackerHitLinkCollection
 *
 * <br><h4>Main input/output collection properties</h4>
 * @param InputClusterCollection Input gas-cluster collection,
 * type edm4hep::SimTrackerHitCollection. <br>
 *
 * @param OutputDigiHitCollection Output direct digitized hit collection,
 * type edm4hep::SenseWireHitCollection. <br>
 *
 * @param DigiSimAssociationCollection Output truth-link collection between
 * final SenseWireHits and contributing gas-cluster SimTrackerHits,
 * type edm4hep::TrackerHitSimTrackerHitLinkCollection. <br>
 *
 * <br><h4>Main configurable properties</h4>
 * @param DetectorName Detector/subdetector name. <br>
 *
 * @param zResolution_mm Gaussian smearing along the sense wire, in mm. <br>
 *
 * @param xyResolution_mm Gaussian smearing perpendicular to the sense wire, in mm. <br>
 *
 * @param XTFileName ROOT file containing the x-t mean/sigma lookup surfaces
 * used by DCHXT2DLUT. <br>
 *
 * @param create_debug_histograms Optional flag to enable QA/debug histograms
 * registered to THistSvc. Output ROOT file name is configured from the Gaudi
 * steering, for example in runDCHdigiV3.py. <br>
 *
 * @param GeoSvcName Geometry service name. <br>
 *
 * @param uidSvcName UniqueIDGenSvc instance name used for reproducible seeds. <br>
 */

#pragma once

// Gaudi Transformer baseclass headers
#include "Gaudi/Property.h"
#include "Gaudi/Accumulators.h"
#include "k4FWCore/Transformer.h"
#include <GaudiKernel/SmartIF.h>
#include "GaudiKernel/ITHistSvc.h"

// Gaudi services
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

// EDM4HEP
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "edm4hep/utils/vector_utils.h"

// EDM4hep wire-hit data model
#include "edm4hep/SenseWireHitCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"

// DD4hep
#include "DDSegmentation/BitFieldCoder.h"
// DD4hep detector extension
#include "detectorCommon/WireTracker_info.h"

// Drift time lookup table
#include "DCHXT2DLUT.h"

// ROOT headers
#include "TH1F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TRandom3.h"
#include <Math/Vector3D.h>

// STL
#include <memory>
#include <string>
#include <vector>

struct WireTrackerDigi_v02 final
    : k4FWCore::MultiTransformer<
          std::tuple<
	  edm4hep::SenseWireHitCollection,
	  edm4hep::TrackerHitSimTrackerHitLinkCollection>(
	  const edm4hep::SimTrackerHitCollection&,
	  const edm4hep::EventHeaderCollection&)> {

  WireTrackerDigi_v02(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  StatusCode finalize() override;

  std::tuple<
	  edm4hep::SenseWireHitCollection,
	  edm4hep::TrackerHitSimTrackerHitLinkCollection>
  operator()(const edm4hep::SimTrackerHitCollection&,
		  const edm4hep::EventHeaderCollection&) const override;

private:
  /// conversion factor mm to cm, static to the class to avoid clash with DD4hep
  static constexpr double MM_TO_CM = 0.1;

  //------------------------------------------------------------------
  //          machinery for geometry

  /// Geometry service name
  Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};

  /// Unique ID service name
  Gaudi::Property<std::string> m_uidSvcName{this, "uidSvcName", "uidSvc", "The name of the UniqueIDGenSvc instance"};

  /// Detector name
  Gaudi::Property<std::string> m_detectorName{this, "DetectorName", "DCH_v2", "Name of the detector"};

  /// Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;

  /// Decoder for the cellID
  const dd4hep::DDSegmentation::BitFieldCoder* m_decoder{nullptr};

  // Geometry helper obtained from the wire-tracker DD4hep extension
  dd4hep::rec::DCH_info* dch_data{nullptr};
  //------------------------------------------------------------------

  // position smearing
  /// Position resolution along the sense wire [mm]
  Gaudi::Property<float> m_z_resolution_mm{
      this, "zResolution_mm", 1.0,
      "Spatial resolution in the z direction (from reading out the wires at both sides) in mm. Default 1 mm."};

  /// Position resolution perpendicular to the sense wire [mm]
  Gaudi::Property<float> m_xy_resolution_mm{this, "xyResolution_mm", 0.1,
                                         "Spatial resolution in the xy direction in mm. Default 0.1 mm."};

  /// Service used to build reproducible event-based random seeds
  SmartIF<IUniqueIDGenSvc> m_uidSvc;

  /// Create an event-reproducible TRandom3 engine from the EventHeader
  TRandom3 CreateRandomEngine(const edm4hep::EventHeaderCollection& headers) const;
  //------------------------------------------------------------------

  /// Print algorithm configuration
  void PrintConfiguration(std::ostream& io) const;

  /// Decode superlayer index from the cellID.
  int CalculateSuperLayerFromCellID(dd4hep::DDSegmentation::CellID id) const {
    return m_decoder->get(id, "superlayer");
  }

  /// Sector is fixed to zero for the current standalone DCH geometry.
  int CalculateSectorFromCellID(dd4hep::DDSegmentation::CellID /*id*/) const {
    return 0;
  }

  /// Convert readout layer fields into the global layer index.
  int CalculateLayerFromCellID(dd4hep::DDSegmentation::CellID id) const {
    return dch_data->CalculateILayerFromCellIDFields(m_decoder->get(id, "layer"), m_decoder->get(id, "superlayer"));
  }

  /// Decode cell index within the layer.
  int CalculateNphiFromCellID(dd4hep::DDSegmentation::CellID id) const { return m_decoder->get(id, "nphi"); }

  /// Convert an EDM4hep vector to a ROOT XYZVector and apply the requested unit scale.
  ROOT::Math::XYZVector Convert_EDM4hepVector_to_XYZVector( const edm4hep::Vector3d& v, double scale) const {
    return {v.x * scale, v.y * scale, v.z * scale};
  }

  /// Convert a ROOT XYZVector back to an EDM4hep vector and apply the requested unit scale.
  edm4hep::Vector3d Convert_XYZVector_to_EDM4hepVector(const ROOT::Math::XYZVector& v, double scale) const {
    return {v.x() * scale, v.y() * scale, v.z() * scale};
  }

  //// ROOT file containing the x-t relation
  Gaudi::Property<std::string> m_xtFileName{this, "XTFileName", "par.root",
	  "ROOT file with x-t relation used to convert radius to drift time" };

  //pointer to the LUT
  std::unique_ptr<DCHXT2DLUT> m_xt2d;

  // Enable diagnostic histograms registered to THistSvc
  Gaudi::Property<bool> m_create_debug_histos{this, "create_debug_histograms", false,
                                              "Enable diagnostic histograms via THistSvc"};

  /// histogram to store distance from sim hit position to the sense wire
  SmartIF<ITHistSvc> m_histSvc;

  /// Cluster electron multiplicity
  TH1F* hNe{nullptr};
  /// Sampled drift-time distribution [ns]
  TH1F* hTd{nullptr};
  /// Drift time versus transverse distance to wire
  TProfile* hXT{nullptr};
  /// Drift time map in the local transverse wire plane
  TProfile2D* hXYT{nullptr};

};

DECLARE_COMPONENT(WireTrackerDigi_v02);


