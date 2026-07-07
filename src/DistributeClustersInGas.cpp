#include "DistributeClustersInGas.h"

// Internal utility headers
#include "BetheBloch.h"
#include "DCHUtilities.h"

// DD4hep and related geometry headers
#include "DD4hep/Detector.h"
#include "DDSegmentation/BitFieldCoder.h"
#include "DD4hep/DetElement.h"

// EDM4hep data model headers
#include <GaudiKernel/MsgStream.h>
#include "GaudiKernel/ITHistSvc.h"

// ROOT headers
#include "TRandom3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

// STL
#include <Math/Vector3D.h>
#include <unordered_map>
#include <cmath>
#include <algorithm>

////////////////////////////////////////////////////////////////////////////
// 							Constructor 							      //
////////////////////////////////////////////////////////////////////////////

// Configure input Geant4 steps and output gas cluster collection.
DistributeClustersInGas::DistributeClustersInGas(
		const std::string& name,ISvcLocator* svcLoc): MultiTransformer(name, svcLoc,
			{
			// Input SimTrackerHit collection and header
			KeyValues("InputSimHitCollection", {"DCHCollection"}),
			KeyValues("HeaderName", {"EventHeader"}),
			},
			{
			// Output cluster-level SimTrackerHits: one object per generated gas cluster
			KeyValues("OutputClusterCollection", {"GasIonizationClusters"}),
      			}) 
{
	// Access DD4hep geometry service.
	m_geoSvc = serviceLocator()->service(m_geoSvcName);
	// Access UniqueIDGenSvc for event-reproducible random seeding
	m_uidSvc = serviceLocator()->service(m_uidSvcName);		
}

//////////////////////////////////////////////////////////////////////////////
//								Initialize Block						    //
//////////////////////////////////////////////////////////////////////////////
StatusCode DistributeClustersInGas::initialize() {
	if (!m_uidSvc) {
		error() << "Unable to get UniqueIDGenSvc" << endmsg;
		return StatusCode::FAILURE;
	}

	if (!m_geoSvc) {
		error() << "Unable to get GeoSvc" << endmsg;
		return StatusCode::FAILURE;
	}

	// Use configured detector name to retrieve the DD4hep detector element.
	const std::string detectorName(m_detectorName.value());
	if (m_geoSvc->getDetector()->detectors().count(detectorName) == 0) {
		error() << "Detector <<" << detectorName << ">> does not exist." << endmsg;
		return StatusCode::FAILURE;
	}

	// Access the DD4hep detector element
	dd4hep::DetElement detectorElement = m_geoSvc->getDetector()->detectors().at(detectorName);

	// Retrieve the wire-tracker geometry extension.
	auto* wireTracker_info = detectorElement.extension<dd4hep::rec::WireTracker_info_struct>();
	if (!wireTracker_info) {
		error() << "No WireTracker_info_struct extension was attached to detector <<"
			<< detectorName << ">>." << endmsg;
		return StatusCode::FAILURE;
	}

	// Use Detector specific helper methods
	// In this case Drift Chamber Helper methods
	dch_data = dynamic_cast<dd4hep::rec::DCH_info*>(wireTracker_info);
	if (!dch_data || !dch_data->IsValid()) {
		error() << "No valid wire-tracker geometry extension was found for detector <<"
				<< detectorName << ">>." << endmsg;
		return StatusCode::FAILURE;
	}

	// Retrieve the sensitive detector associated with DCH
	dd4hep::SensitiveDetector sensitiveDetector = m_geoSvc->getDetector()->sensitiveDetector(detectorName);
	if (!sensitiveDetector.isValid()) {
		error() << "No valid Sensitive Detector was found for detector <<" 
			<< detectorName << ">>." << endmsg;
		return StatusCode::FAILURE;
	}

	// Access the detector readout
	dd4hep::Readout readout = sensitiveDetector.readout();

	// Retrieve the CellID decoder from the readout
	m_decoder = readout.idSpec().decoder();
	if (!m_decoder) {
		error() << "Failed to retrieve CellID decoder from detector readout." << endmsg;
		return StatusCode::FAILURE;
	}

	info() << name() << " initialized for detector " << detectorName << endmsg;

	// Reset the internal event counter at job start
	m_event_counter.reset();
//---------------------------------------------------------------------------------------------//

// Create and register optional debug histograms
	if (m_create_debug_histos.value()) {
		m_histSvc = service("THistSvc");
		if (!m_histSvc) {
			error() << "Unable to locate Histogram Service" << endmsg;
			return StatusCode::FAILURE;
		}

		// Histograms for generated cluster positions, and energies
		hX = new TH1F("hX",
				"Generated cluster X positions;X [mm];Counts",
				50, 0.0, 2500.0);
		if (m_histSvc->regHist("/rec/hX_clusters", hX).isFailure()) {
			error() << "Couldn't register hX histogram" << endmsg;
			return StatusCode::FAILURE;
		}
		hY = new TH1F("hY",
                                "Generated cluster Y positions;Y [mm];Counts",
                                50, 0.0, 2500.0);
                if (m_histSvc->regHist("/rec/hY_clusters", hY).isFailure()) {
                        error() << "Couldn't register hY histogram" << endmsg;
                        return StatusCode::FAILURE;
                }
		hZ = new TH1F("hZ",
                                "Generated cluster Z positions;Z [mm];Counts",
                                50, 0.0, 2500.0);
                if (m_histSvc->regHist("/rec/hZ_clusters", hZ).isFailure()) {
                        error() << "Couldn't register hZ histogram" << endmsg;
                        return StatusCode::FAILURE;
                }
		hXYZ = new TH3F("hXYZ",
				"Generated cluster positions;X [mm];Y [mm];Z [mm]",
				50, 0.0, 2500.0,
				50, 0.0, 2500.0,
				50, 0.0, 2500.0);
		if (m_histSvc->regHist("/rec/hXYZ_clusters", hXYZ).isFailure()) {
			error() << "Couldn't register hXYZ histogram" << endmsg;
			return StatusCode::FAILURE;
		}
		heDep = new TH1F("heDep",
				"Assigned cluster energy;E_{cluster} [GeV];Counts",
				50, 0.0, 1.0e-5);
		if (m_histSvc->regHist("/rec/heDep_clusters", heDep).isFailure()) {
			error() << "Couldn't register heDep histogram" << endmsg;
			return StatusCode::FAILURE;
		}
		hR = new TH1F("hR",
				"Generated cluster radius;R [mm];Counts",
				50, 0.0, 2500.0);
		if (m_histSvc->regHist("/rec/hR_clusters", hR).isFailure()) {
			error() << "Couldn't register hR histogram" << endmsg;
			return StatusCode::FAILURE;
		}
		hRvsZ = new TH2F("hRvsZ",
				"Generated cluster R vs Z;Z [mm];R [mm]", 100, 0.0, 2500.0, 100, 0.0, 2500.0);
				if (m_histSvc->regHist("/rec/hRvsZ_clusters", hRvsZ).isFailure()) {
					error() << "Couldn't register hRvsZ histogram" << endmsg;
					return StatusCode::FAILURE;
				}
		hPathLength = new TH1F("hPathLength", "Path Length per Geant4 step;Path Length [mm];Entries",
			50, 0.0, 12.0);
			if (m_histSvc->regHist("/rec/hPathLength_clusters", hPathLength).isFailure()) {
				error() << "Couldn't register hPathLength histogram" << endmsg;
				return StatusCode::FAILURE;
			}
		hPLvsnC = new TH2F("hPLvsnC",
			"Path Length vs nCell;PathLength [mm];nCell",
			50, 0.0, 12.0,
            50, 0.0, 150.0);
			if (m_histSvc->regHist("/rec/hPLvsnC_clusters", hPLvsnC).isFailure()) {
				error() << "Couldn't register hPLvsnC histogram" << endmsg;
                return StatusCode::FAILURE;
            }

		// Histograms for cluster and inter-cluster spacing	
        hNcl = new TH1F("hNcl", "n_{cl} per cm;n_{cl} / cm;Entries", 25, 0.0, 50.0);
		if (m_histSvc->regHist("/rec/hNcl_clusters", hNcl).isFailure()) {
            error() << "Couldn't register hNcl histogram" << endmsg;
            return StatusCode::FAILURE;
        }
		hClSpacing_mm = new TH1F("hClSpacing_mm", "Inter-cluster spacing;spacing [mm];Entries",
			50, 0.0, 5.0);
			if (m_histSvc->regHist("/rec/hClSpacing_mm_clusters", hClSpacing_mm).isFailure()) {
                error() << "Couldn't register hClSpacing_mm histogram" << endmsg;
                return StatusCode::FAILURE;
            }

	} // End of histogram registration

	// Finish successful initialization
	return StatusCode::SUCCESS;

} // End of Initialize Block

// Build an event-reproducible random-number generator from the event header
TRandom3
DistributeClustersInGas::CreateRandomEngine(
		const edm4hep::EventHeaderCollection& headers) const {
      	const auto seed = m_uidSvc->getUniqueID(headers, this->name());
      	TRandom3 rng(seed);
      	return rng;
}

//===============================================================================//
//                        Cluster density from Bethe-Bloch                       //
//===============================================================================//
// Estimate the mean number of primary ionization clusters per cm 
// for a given particle beta*gamma and mass
double DistributeClustersInGas::dNcldx(double betagamma,const edm4hep::MCParticle& mc) const {
	if (betagamma <= 0.0 || !std::isfinite(betagamma)) {
		return 0.0;
	}
	if (!mc.isAvailable()) {
		return 0.0;
	}
	// Get particle momentum from beta*gamma and MC mass
	const double m_GeV = mc.getMass();
	if (m_GeV <= 0.0 || !std::isfinite(m_GeV)) {
		return 0.0;
	}
	const double p_GeV = betagamma * m_GeV;
	const double m_MeV = m_GeV * 1.0e3;
	const double I_MeV = m_MeanExcEnergy_eV.value() * 1.0e-6;

	// Evaluate Bethe-Block
	const double dEdx_MeVcm2_per_g = BB::bethe_bloch(p_GeV, m_MeV, I_MeV);
			
	// Convert to MeV/cm
	const double dEdx_MeV_per_cm = dEdx_MeVcm2_per_g * m_GasDensity_g_cm3.value();

	// Convert deposited energy into mean number of primary ionization
	const double lambda_per_cm = (dEdx_MeV_per_cm * 1.0e6) / m_W_eff_eV.value();

	return (lambda_per_cm > 0.0) ? lambda_per_cm : 0.0;
}

// look up the MCParticle production vertex and 
// check if it is located inside the active gas volume
bool DistributeClustersInGas::isParticleCreatedInsideActiveGasVolume(
	const edm4hep::MCParticle& mc) const {
		if (!dch_data || !mc.isAvailable())  return false;

		// Use MC production vertex to decide whether a secondary was born in active gas.
		const auto vertex = mc.getVertex();  // in mm

		// Compute squared radial distance of the vertex from the detector axis
		const double vertexRsquared = vertex[0] * vertex[0] + vertex[1] * vertex[1];
		const double vertexZabs = std::fabs(vertex[2]);

		// Read the active gas volume half length, inner radius^2, and outer radius^2 in mm
		const double detectorHalfLength_mm = dch_data->Lhalf / dd4hep::mm;
		const double innerRadiusSquared_mm2 = std::pow(dch_data->rin / dd4hep::mm, 2);
		const double outerRadiusSquared_mm2 = std::pow(dch_data->rout / dd4hep::mm, 2);

		// Compare production point with active gas z and radial boundaries.
		return (vertexZabs < detectorHalfLength_mm) &&
			(vertexRsquared > innerRadiusSquared_mm2) &&
			(vertexRsquared < outerRadiusSquared_mm2);
}

////////////////////////////////////////////////////////////////////////////////
//                            Operator Block                                  //
////////////////////////////////////////////////////////////////////////////////
std::tuple<edm4hep::SimTrackerHitCollection>
DistributeClustersInGas::operator()(const edm4hep::SimTrackerHitCollection& input_sim_hits,
                                    const edm4hep::EventHeaderCollection& headers) const
{
	if (m_event_counter.value() % 10 == 0) {
		info() << "Processing Events: " << m_event_counter.value() << endmsg;
       	
	}
	
	//Reproducible event-level random engine
	auto rng = CreateRandomEngine(headers);

	// Output collection:
	// one output object per generated gas cluster
	edm4hep::SimTrackerHitCollection output_cluster_hits;
	
	// Keep accepted Geant4 steps and merge them at cell level
    struct StepSegment {
	    edm4hep::SimTrackerHit simHit;     // original Geant4 step
	    edm4hep::MCParticle mcParticle;    // particle producing this step

	    ROOT::Math::XYZVector stepMid_cm;               // midpoint of the Geant4 step [cm]
	    ROOT::Math::XYZVector trackDir;                 // unit direction of the Geant4 step
	    double stepLength_mm = 0.0;        // path length of this step [mm]

	    double time_ns = 0.0;              // Geant4 hit time
	    double eDep_GeV = 0.0;             // deposited energy in this step
    };

	// Merge all accepted Geant4 steps belonging to the same readout cell.
    struct CellAccum {
	    std::vector<StepSegment> segments; // all accepted Geant4 steps in this cell
	    double totalPathLength_mm = 0.0;   // merged path length in this cell
	    double totalEDep_GeV = 0.0;        // summed deposited energy in this cell
    };

    std::unordered_map<uint64_t, CellAccum> cellMap;

    // loop over all Geant4 steps in the input collection
    int loop_index = 0;
    for (const auto& input_sim_hit : input_sim_hits) 
	{
	    ++loop_index;

	    const auto mcParticle = input_sim_hit.getParticle();
	    const int ilayer = calculateLayerFromCellID(input_sim_hit.getCellID());
	    const int nphi   = calculateNphiFromCellID(input_sim_hit.getCellID());
		const uint64_t cellID = static_cast<uint64_t>(input_sim_hit.getCellID());

	    // Hit selection criteria:
	    const bool isPrimaryHit = mcParticle.isAvailable() &&
		mcParticle.getParents().empty() &&
		!input_sim_hit.isProducedBySecondary();

		// Secondary step accepted only if created inside the active gas volume
	    const bool isSecondaryInsideActiveGasVolume =
		mcParticle.isAvailable() &&
		!mcParticle.getParents().empty() &&
		input_sim_hit.isProducedBySecondary() &&
		isParticleCreatedInsideActiveGasVolume(mcParticle);

	    if (m_keepOnlyPrimaryHits.value()) {
			if (!isPrimaryHit) continue;
		} 
		else {
			if (!(isPrimaryHit || isSecondaryInsideActiveGasVolume)) continue;
		}
    
		// Reject null or unphysical Geant4 step lengths
		const double stepLength_mm = input_sim_hit.getPathLength();
		//if (stepLength_mm <= 0.0) continue;
    
		// Build normalized Geant4 step direction from momentum
	    ROOT::Math::XYZVector trackDir(
			input_sim_hit.getMomentum().x,
		    input_sim_hit.getMomentum().y,                     
			input_sim_hit.getMomentum().z
		);
		
		if (std::sqrt(trackDir.Dot(trackDir)) <= 0.0) continue;
	    trackDir = trackDir.Unit();

		// Convert Geant4 step position from mm to cm
	    const auto pos = input_sim_hit.getPosition();
	    ROOT::Math::XYZVector stepMid_cm(
			pos.x * MM_TO_CM,
		    pos.y * MM_TO_CM,
			pos.z * MM_TO_CM
		);

	    if (m_create_debug_histos.value()) {
			hPathLength->Fill(stepLength_mm);
		  	hPLvsnC->Fill(stepLength_mm, nphi);
	    }

	    // Append this accepted Geant4 step to the corresponding cell
	    auto& cellAccum = cellMap[cellID];

		// Cache the accepted Geant4 step information needed for per-cell cluster generation.
	    StepSegment seg;
	    seg.simHit        = input_sim_hit;
		seg.mcParticle    = mcParticle;
	    seg.stepMid_cm    = stepMid_cm;
	    seg.trackDir      = trackDir;
	    seg.stepLength_mm = stepLength_mm;
	    seg.time_ns       = input_sim_hit.getTime();
	    seg.eDep_GeV      = input_sim_hit.getEDep();

	    cellAccum.segments.push_back(seg);
	    cellAccum.totalPathLength_mm += stepLength_mm;
	    cellAccum.totalEDep_GeV      += input_sim_hit.getEDep();

	    debug () << "Accepted step #" << loop_index 
				<< " layer=" << ilayer
	    		<< " cell=" << nphi
	    		<< " pathLength_mm=" << stepLength_mm
				<< " TotalPathLength_mm=" << cellAccum.totalPathLength_mm
	    		//<< " eDep_GeV=" << input_sim_hit.getEDep()
				<< endmsg;
	
	} // End of loop over simulated Hits

	//===============================================================================//
	//                           Generate Cluster per Cell                           //
	//===============================================================================//
	// mean number of electrons per cluster for the He-based gas mixture
	static const double Ne_mean = mean_cluster_size_He();
	
	// loop over cellMap to generate clusters for each cell
	for (const auto& kv : cellMap) 
	{
	    const uint64_t cellID = kv.first;
	    const CellAccum& cell = kv.second;

	    if (cell.segments.empty()) continue;
	    if (cell.totalPathLength_mm <= 0.0) continue;

		// Convert total path length from mm to cm
	    const double totalL_mm = cell.totalPathLength_mm;
	    const double totalL_cm = totalL_mm * MM_TO_CM;

	    // Compute mean number of clusters
	    double mu_clusters = 0.0;
		//loop over each segment
	    for (const auto& seg : cell.segments) {
			if (seg.stepLength_mm <= 0.0) continue;
			if (!seg.mcParticle.isAvailable()) continue;

			// Particle beta*gamma used in Bethe-Bloch cluster-density estimate
			const double bg = compute_beta_gamma(seg.mcParticle);
			const double lambda_electrons_per_cm = dNcldx(bg, seg.mcParticle);

			// Convert mean primary electrons/cm to mean physical clusters/cm
			const double lambda_clusters_per_cm = lambda_electrons_per_cm / Ne_mean;
			mu_clusters += lambda_clusters_per_cm * (seg.stepLength_mm * MM_TO_CM);
	    }

		// Distribute Cluster using Trancated Poisson Distribution
	    const int nClusters = ZT_poisson(mu_clusters, rng);
	    if (nClusters <= 0) continue;

	    if (m_create_debug_histos.value()) {
			const double nClusters_per_cm = (totalL_cm > 0.0) ? 
			static_cast<double>(nClusters) / totalL_cm : 0.0;
			hNcl->Fill(nClusters_per_cm);
	    }
		//------------------------------------------------------------------------//

	    // Build cumulative step-end positions along the merged cell path.
	    std::vector<double> stepEndPositions_mm;
	    stepEndPositions_mm.reserve(cell.segments.size());

		// Running end position along the merged path of this cell [mm].
	    double pathEnd_mm = 0.0;
		// Convert step lengths into cumulative boundaries for later cluster-to-step mapping.
	    for (const auto& segment : cell.segments) {
			// Move the path end by the length of the current Geant4 step.
		  	pathEnd_mm += segment.stepLength_mm;
			// Store where this Geant4 step ends along the merged cell path.
		  	stepEndPositions_mm.push_back(pathEnd_mm);
	    }
		//-----------------------------------------------------------------------//

		// Generate ordered cluster positions along the merged path of this cell
		const std::vector<double> clusterPositions_mm = generate_cluster_positions_mm(
			totalL_mm, nClusters, rng);

		// Fill inter-cluster spacing from consecutive ordered positions
		if (m_create_debug_histos.value() && hClSpacing_mm) {
			for (size_t i = 1; i < clusterPositions_mm.size(); ++i) {
				const double spacing_mm = clusterPositions_mm[i] - clusterPositions_mm[i - 1];
				if (spacing_mm >= 0.0) {
					hClSpacing_mm->Fill(spacing_mm);
				}
			}
		}
		//-------------------------------------------------------------------------//

		// Create one EDM cluster object per generated cluster
		for (size_t icl = 0; icl < clusterPositions_mm.size(); ++icl) {

		  	// Global cluster coordinate along the merged cell path [mm]
			const double sGlobal_mm = clusterPositions_mm[icl];

			// Find the first Geant4 step whose end position is after this cluster position.
		  	auto it = std::lower_bound(
				stepEndPositions_mm.begin(),
			 	stepEndPositions_mm.end(), 
				sGlobal_mm
			);

			// Skip this cluster if no matching Geant4 step is found.
		  	if (it == stepEndPositions_mm.end()) {
				warning() << "Generated cluster position is outside the accumulated step path: "
				<< "sGlobal_mm = "<< sGlobal_mm
				<< ", lastStepEnd_mm=" << stepEndPositions_mm.back()
				<< ", cellID=" << cellID
				<< endmsg;
				continue;
			}

			// Index of the Geant4 step that contains this generated cluster.
		  	const size_t iseg = std::distance(stepEndPositions_mm.begin(), it);
			// Original Geant4 step used to reconstruct the 3D cluster position.
		  	const auto& seg = cell.segments[iseg];

			// End position of the previous Geant4 step along the merged path [mm].
		  	const double previousStepEnd_mm = (iseg == 0) ? 0.0 : stepEndPositions_mm[iseg - 1];
			// Local position of this cluster inside the selected Geant4 step [mm].
		  	const double sLocal_mm  = sGlobal_mm - previousStepEnd_mm;
			// Convert local cluster position from mm to cm.
		  	const double sLocal_cm  = sLocal_mm * MM_TO_CM;

		  	// Reconstruct start point of the selected Geant4 step [cm].
		  	const ROOT::Math::XYZVector stepStart_cm =
		      		seg.stepMid_cm - 0.5 * (seg.stepLength_mm * MM_TO_CM) * seg.trackDir;

		  	// Place the cluster along the selected Geant4 step
		  	const ROOT::Math::XYZVector clusterPos_cm = stepStart_cm + sLocal_cm * seg.trackDir;

			// Store cluster position in EDM4hep units [mm].
		  	edm4hep::Vector3d clusterPos_mm {
				clusterPos_cm.X() / MM_TO_CM,
			    clusterPos_cm.Y() / MM_TO_CM,
			    clusterPos_cm.Z() / MM_TO_CM
		  	};
			//-------------------------------------------------------------------------//

		  	// Share total cell energy uniformly over generated clusters
		  	const float clusterEDep = static_cast<float>(
				cell.totalEDep_GeV / static_cast<double>(nClusters)
			);

		  	if (m_create_debug_histos.value()) {
				const double r_mm = std::hypot(clusterPos_mm.x, clusterPos_mm.y);
				hX->Fill(clusterPos_mm.x);
				hY->Fill(clusterPos_mm.y);
				hZ->Fill(clusterPos_mm.z);
				hXYZ->Fill(clusterPos_mm.x, clusterPos_mm.y, clusterPos_mm.z);
				heDep->Fill(clusterEDep);
				hR->Fill(r_mm);
				hRvsZ->Fill(clusterPos_mm.z, r_mm);
		  	}


			//===============================================================================//
			//          Store the generated gas cluster in the EDM output collection         //
			//===============================================================================//
		  	auto out_hit = output_cluster_hits.create();
		  	out_hit.setCellID(cellID);
		  	out_hit.setPosition(clusterPos_mm);
		  	out_hit.setTime(static_cast<float>(seg.time_ns));
		  	out_hit.setEDep(clusterEDep);

			// Preserve the source MCParticle relation through the standard SimTrackerHit field.
			out_hit.setParticle(seg.mcParticle);

			// Keep original step direction information where possible.
			out_hit.setMomentum(seg.simHit.getMomentum());
			//out_hit.setPathLength(0.0f);

		  	debug() << "Cell " << cellID
		  		<< " cluster #" << icl
		  		<< " segment=" << iseg
		  		<< " sLocal_mm=" << sLocal_mm
		  		<< " pos_mm=("
		  		<< clusterPos_mm.x << ", "
		  		<< clusterPos_mm.y << ", "
		  		<< clusterPos_mm.z << ")"
		  		<< endmsg;
	    	
		} // End of loop over generated clusters in this cell

    } // End of the loop over Cell

	// Events counter
	++m_event_counter;

	// Return one cluster-level SimTrackerHit per generated gas cluster.
	return std::make_tuple(std::move(output_cluster_hits));

} // End of the operator block
