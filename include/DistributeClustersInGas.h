#pragma once

/** ================= DistributeClustersInGas =================
 * Gaudi Algorithm for gas-ionization cluster generation in a wire tracker.
 *
 * @author Muhammad Saiel
 * @date   2026-04
 *
 * <h4>Purpose</h4>
 * This algorithm is the upstream gas-cluster generation stage of the wire-tracker
 * digitization chain. It converts Geant4 tracker-step information
 * (edm4hep::SimTrackerHitCollection) into a new edm4hep::SimTrackerHitCollection
 * where each output SimTrackerHit represents one generated gas-ionization cluster.
 *
 * This code uses DD4hep natural length units (cm), while EDM4hep positions and
 * path lengths are stored in mm. Please handle unit conversions carefully.
 *
 * <h4>Scope of this stage</h4>
 * This module is responsible for:
 *   - selecting the Geant4 SimTrackerHits used for gas-cluster generation
 *   - estimating the cluster density from the gas model
 *   - sampling the number of clusters produced in each Geant4 step/cell
 *   - sampling cluster positions along the Geant4 step
 *   - storing one output edm4hep::SimTrackerHit per generated gas cluster
 *
 * <h4>Input</h4>
 * @param InputSimHitCollection Input Geant4 tracker-hit collection,
 * type edm4hep::SimTrackerHitCollection.
 *
 * @param HeaderName Event header collection used to build reproducible random
 * seeds.
 *
 * <h4>Output</h4>
 * @param OutputClusterCollection Output gas-cluster collection,
 * type edm4hep::SimTrackerHitCollection.
 *
 * <h4>Main configurable properties</h4>
 * @param DetectorName Detector/subdetector name in DD4hep.
 *
 * @param MeanExcitationEnergy_eV Mean excitation energy of the gas mixture.
 *
 * @param GasDensity_g_cm3 Gas density used to convert Bethe-Bloch stopping power
 * into energy loss per unit length.
 *
 * @param W_eff_eV Effective energy per ionization cluster used to convert energy
 * loss into mean cluster density.
 *
 * @param KeepOnlyPrimaryHits If true, only primary Geant4 hits are converted into
 * gas clusters. If false, primary hits and selected secondary hits created inside
 * the active gas volume are clusterized.
 *
 * @param create_debug_histograms Optional flag to enable diagnostic histograms
 * via THistSvc.
 */

#include "Gaudi/Accumulators.h"
#include "Gaudi/Property.h"
#include "k4FWCore/Transformer.h"

#include <GaudiKernel/SmartIF.h>
#include "GaudiKernel/ITHistSvc.h"

#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticle.h"
#include "edm4hep/SimTrackerHitCollection.h"

#include "DDRec/DCH_info.h"
#include "DDSegmentation/BitFieldCoder.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TRandom3.h"

#include <string>

struct DistributeClustersInGas final
    : k4FWCore::MultiTransformer<
          std::tuple<edm4hep::SimTrackerHitCollection>(
              const edm4hep::SimTrackerHitCollection&,
              const edm4hep::EventHeaderCollection&)> {

  DistributeClustersInGas(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;

  std::tuple<edm4hep::SimTrackerHitCollection>
  operator()(const edm4hep::SimTrackerHitCollection& inputSimHits,
             const edm4hep::EventHeaderCollection& headers) const override;

private:
  /// Conversion factor from EDM4hep length unit (mm) to DD4hep natural unit (cm)
  static constexpr double MM_TO_CM = 0.1;

  // --------------------------------------------------------------------------
  // Geometry and services
  // --------------------------------------------------------------------------

  /// Geometry service name
  Gaudi::Property<std::string> m_geoSvcName{ this,
    "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};

  /// Unique ID service name
  Gaudi::Property<std::string> m_uidSvcName{this, 
    "uidSvcName", "uidSvc", "The name of the UniqueIDGenSvc instance"};

  /// Detector/subdetector name in DD4hep
  Gaudi::Property<std::string> m_detectorName{this, 
    "DetectorName", "DCH_v2", "Name of the detector"};

  /// Pointer to geometry service
  SmartIF<IGeoSvc> m_geoSvc;

  /// Service used to build reproducible event-based random seeds
  SmartIF<IUniqueIDGenSvc> m_uidSvc;

  /// CellID decoder from readout
  const dd4hep::DDSegmentation::BitFieldCoder* m_decoder{nullptr};

  /// Detector geometry extension.
  /// TODO: replace DCH_info by the common WireTracker_info once available.
  dd4hep::rec::DCH_info* dch_data{nullptr};

  /// Event counter
  mutable Gaudi::Accumulators::Counter<
      Gaudi::Accumulators::atomicity::full, unsigned int>
      m_event_counter{this, "EventsProcessed"};

  int calculateLayerFromCellID(dd4hep::DDSegmentation::CellID id) const {
    return dch_data->CalculateILayerFromCellIDFields(
        m_decoder->get(id, "layer"),
        m_decoder->get(id, "superlayer"));
  }

  int calculateNphiFromCellID(dd4hep::DDSegmentation::CellID id) const {
    return m_decoder->get(id, "nphi");
  }

  // --------------------------------------------------------------------------
  // Gas / ionization properties
  // --------------------------------------------------------------------------

  /// Mean excitation energy of the gas
  Gaudi::Property<double> m_MeanExcEnergy_eV{this, 
    "MeanExcitationEnergy_eV", 48.48,"Mean excitation energy I in eV for the gas mixture"};

  /// Gas density used to convert dE/dx into MeV/cm
  Gaudi::Property<double> m_GasDensity_g_cm3{this, 
    "GasDensity_g_cm3", 3.984e-4, "Gas density in g/cm3 used to convert MeV.cm2/g into MeV/cm"};

  /// Effective energy per ionization
  Gaudi::Property<double> m_W_eff_eV{this, 
    "W_eff_eV", 35.0, "Effective energy per ion pair in eV"};

  /// If true, only Geant4 hits from the primary particle are converted into gas clusters.
  Gaudi::Property<bool> m_keepOnlyPrimaryHits{this, 
    "KeepOnlyPrimaryHits", true,
    "If true, keep only primary Geant4 SimTrackerHits for gas-cluster generation"};

  /// Build an event-reproducible random engine from the EventHeader
  TRandom3 CreateRandomEngine(const edm4hep::EventHeaderCollection& headers) const;

  /// Convert Bethe-Bloch + gas properties into mean number of clusters per cm
  double get_dNcldx_per_cm(double betagamma, const edm4hep::MCParticle& mc) const;

  /// Check whether a secondary particle was produced inside the active gas volume
  bool isParticleCreatedInsideActiveGasVolume(const edm4hep::MCParticle& mc) const;

  /// Enable  histograms registered to THistSvc
  Gaudi::Property<bool> m_create_debug_histos{this, 
    "create_debug_histograms", false,"Enable diagnostic histograms for gas-cluster generation via THistSvc"};

  /// Histogram service
  SmartIF<ITHistSvc> m_histSvc;

  /// Geant4 step path length
  TH1F* hPathLength{nullptr};

  /// Path length versus cell index
  TH2F* hPLvsnC{nullptr};

  /// Number of generated clusters per cm
  TH1F* hNcl{nullptr};

  /// Inter-cluster spacing in mm
  TH1F* hClSpacing_mm{nullptr};

  /// Generated cluster X positions [mm]
  TH1F* hX{nullptr};

  /// Generated cluster Y positions [mm]
  TH1F* hY{nullptr};

  /// Generated cluster Z positions [mm]
  TH1F* hZ{nullptr};

  /// Generated cluster 3D positions [mm]
  TH3F* hXYZ{nullptr};

  /// Energy assigned to each generated cluster [GeV]
  TH1F* heDep{nullptr};

  /// Generated cluster radial distance from detector axis [mm]
  TH1F* hR{nullptr};

  /// Generated cluster radial distance versus z [mm]
  TH2F* hRvsZ{nullptr};
};

DECLARE_COMPONENT(DistributeClustersInGas)