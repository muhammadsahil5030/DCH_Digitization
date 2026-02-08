/** ======= DCHdigi_v01 ==========
 * Gaudi Algorithm for DCH digitization
 *
 *
 * @author Alvaro Tolosa-Delgado, Brieuc Francois
 * @date   2024-08
 *
 * <h4>Input collections and prerequisites</h4>
 * Processor requires a collection of SimTrackerHits <br>
 * This code uses DD4hep length natural unit (cm), but EDM4hep data is (usually) in mm. Please be careful with units.
 * <br> <h4>Output</h4> Processor produces collection of Digitized hits of Drift Chamber v2<br>
 * @param DCH_simhits The name of input collection, type edm4hep::SimTrackerHitCollection <br>
 * (default name empty) <br>
 * @param DCH_DigiCollection The name of out collection, type extension::SenseWireHitCollection <br>
 * (default name DCH_DigiCollection) <br>
 * @param DCH_name DCH subdetector name <br>
 * (default value DCH_v2) <br>
 * @param calculate_dndx Optional flag to calcualte dNdx information <br>
 * (default value false) <br>
 * @param fileDataAlg File needed for calculating cluster count and size <br>
 * (default value /eos/.../DataAlgFORGEANT.root) <br>
 * @param zResolution_mm Resolution (sigma for gaussian smearing) along the sense wire, in mm <br>
 * (default value 1 mm) <br>
 * @param xyResolution_mm Resolution (sigma for gaussian smearing) perpendicular the sense wire, in mm <br>
 * (default value 0.1 mm) <br>
 * @param create_debug_histograms Optional flag to create debug histograms <br>
 * (default value false) <br>
 * @param GeoSvcName Geometry service name <br>
 * (default value GeoSvc) <br>
 * @param uidSvcName The name of the UniqueIDGenSvc instance, used to create seed for each event/run, ensuring
 * reproducibility. <br> (default value uidSvc) <br> <br>
 */

#ifndef DCHDIGI_V01_H
#define DCHDIGI_V01_H

// Gaudi Transformer baseclass headers
#include "Gaudi/Property.h"
#include "k4FWCore/Transformer.h"

// Gaudi services
#include "k4Interface/IGeoSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

// EDM4HEP
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"

// EDM4HEP extension
#include "extension/SenseWireHitCollection.h"
#include "extension/SenseWireHitSimTrackerHitLinkCollection.h"

// DD4hep
#include "DD4hep/Detector.h" // for dd4hep::VolumeManager
#include "DDSegmentation/BitFieldCoder.h"

// STL
#include <random>
#include <string>

// data extension for detector DCH_v2
#include "DDRec/DCH_info.h"

// ROOT headers
#include "TFile.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TTree.h"
#include "TProfile.h"
#include "TGraph.h"
#include "DCHXT2DLUT.h"
#include <memory>


//requirments for the xtrelTIME
class XTRELTIME;

class TF1;
class TFile;
class TGraph2D;

/// constant to convert from mm (EDM4hep) to DD4hep (cm)

struct DCHdigi_v01 final
    : k4FWCore::MultiTransformer<
          std::tuple<extension::SenseWireHitCollection, extension::SenseWireHitSimTrackerHitLinkCollection>(
              const edm4hep::SimTrackerHitCollection&, const edm4hep::EventHeaderCollection&)> {
  DCHdigi_v01(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  StatusCode finalize() override;

  std::tuple<extension::SenseWireHitCollection, extension::SenseWireHitSimTrackerHitLinkCollection>
  operator()(const edm4hep::SimTrackerHitCollection&, const edm4hep::EventHeaderCollection&) const override;

private:
  /// conversion factor mm to cm, static to the class to avoid clash with DD4hep
  static constexpr double MM_TO_CM = 0.1;

  //------------------------------------------------------------------
  //          machinery for geometry

  /// Geometry service name
  Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};
  Gaudi::Property<std::string> m_uidSvcName{this, "uidSvcName", "uidSvc", "The name of the UniqueIDGenSvc instance"};

  /// Detector name
  Gaudi::Property<std::string> m_DCH_name{this, "DCH_name", "DCH_v2", "Name of the Drift Chamber detector"};

  /// Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;

  /// Decoder for the cellID
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;

  /// Pointer to drift chamber data extension
  dd4hep::rec::DCH_info* dch_data = {nullptr};

  double m_halfChamberLength_mm = 0.0;	// half length along z

  //------------------------------------------------------------------
  //          machinery for smearing the position

  /// along the sense wire position resolution in mm
  Gaudi::Property<float> m_z_resolution{
      this, "zResolution_mm", 1.0,
      "Spatial resolution in the z direction (from reading out the wires at both sides) in mm. Default 1 mm."};
  /// xy resolution in mm
  Gaudi::Property<float> m_xy_resolution{this, "xyResolution_mm", 0.1,
                                         "Spatial resolution in the xy direction in mm. Default 0.1 mm."};

  /// create seed using the uid
  SmartIF<IUniqueIDGenSvc> m_uidSvc;

  /// Create random engine, initialized with seed out of Event Header
  std::tuple<std::mt19937_64, TRandom3> CreateRandomEngines(const edm4hep::EventHeaderCollection& headers) const;
  //------------------------------------------------------------------
  //        ancillary functions

  bool IsFileGood(std::string& ifilename) const { return std::ifstream(ifilename).good(); }

  /// Print algorithm configuration
  void PrintConfiguration(std::ostream& io);

  /// Send error message to logger and then throw exception
  void ThrowException(std::string s) const;

  int CalculateLayerFromCellID(dd4hep::DDSegmentation::CellID id) const {
    // return m_decoder->get(id, "layer") + dch_data->nlayersPerSuperlayer * m_decoder->get(id, "superlayer") + 1;
    return dch_data->CalculateILayerFromCellIDFields(m_decoder->get(id, "layer"), m_decoder->get(id, "superlayer"));
  }

  int CalculateNphiFromCellID(dd4hep::DDSegmentation::CellID id) const { return m_decoder->get(id, "nphi"); }

  TVector3 Convert_EDM4hepVector_to_TVector3(const edm4hep::Vector3d& v, double scale) const {
    return {v[0] * scale, v[1] * scale, v[2] * scale};
  };
  edm4hep::Vector3d Convert_TVector3_to_EDM4hepVector(const TVector3& v, double scale) const {
    return {v.x() * scale, v.y() * scale, v.z() * scale};
  };

  bool IsParticleCreatedInsideDriftChamber(const edm4hep::MCParticle&) const;
  
  //---------------------------Changes made by Muhammad Saiel------------------------------//
  // Ionization due to MC Particle
  Gaudi::Property<double> m_MCParticle{this, "MCParticle", 13,
	  "ID of the monte carlo particle that produce the hit (defaul: Muon)"};
  
  //Mean excitation energy
  Gaudi::Property<double> m_MeanExcEnergy_eV{this, "MeanExcitationEnergy_eV", 48.48, 
	  "Mean excitation Energy I in eV for the gas (Ref: materials_o1_v02.xml)"};

  //Gas density(default: He-isobutane mixture in 90:10)
  Gaudi::Property<double> m_GasDensity_g_cm3{this, "GasDensity_g_cm3", 3.984e-4,
	  "Gas density in g/cm3 used to convert MeV.cm2/g -> MeV/cm"};

  //Energy to produce single ion pair (default: He-isobutane mixture in 90:10)
  Gaudi::Property<double> m_W_eff_eV{this, "W_eff_eV", 35.0,
	  "Effective energy per cluster in eV (W_eff)"};

  //Mass of the MC particle
  Gaudi::Property<double> m_MassForBB_GeV{this, "MassForBB_GeV", 0.105658,
	  "Reference particle mass (GeV) for BetaGamma conversion (default muon mass)"};

  //function to calculate cluster density using Betagamma
  double get_dNcldx_per_cm(double betagamma) const;

  //this block is used for drift time parameterization:
  Gaudi::Property<std::string> m_xtFileName{this, "XTFileName", "par.root",
	  "ROOT file with x-t relation used to convert radius to drift time" };

  mutable XTRELTIME* m_xtHelper{nullptr};
  double electronDriftTime(double r_cm, TRandom3& myRandom) const;
  double electronDriftVelocity_cm_per_us(double r_cm) const;
  double electricField_V_per_cm(double r_cm) const;

  std::unique_ptr<DCHXT2DLUT> m_xt2d;

  //Gaudi properties for the Gas Gain/ polya parameters:
  Gaudi::Property<double> m_GasGain{this, "GasGain", 2.0e5,
	  "Mean gas gain used in Polya distribution (arbitrary units)" };
  
  //Polya distribution shape parameter
  Gaudi::Property<double> m_PolyaTheta{this, "PolyaTheta", 0.5,
	  "Polya shape parameter theta (typicaly around 5 for drift tube)"};
  
  TF1* m_polya{nullptr};

  //Function to produce avalunche charge sample from polya
  double avalancheCharge(TRandom3& myRandom) const;

  //----------------- this block is looking for pulse shaping  --------------------//

  //Scale parameter
  Gaudi::Property<double> m_pulseAmplitudeScale{this, "PulseAmplitudeScale", 5.e-6,
          "Global Scale factor converting avalache charge to voltage"};
  
  //Pulse rise time
  Gaudi::Property<double> m_pulseRiseTime_ns{this, "PulseRiseTime_ns", 1.0,
          "Rise time of the single electron pulse (ns)"};

  //Pulse decay time
  Gaudi::Property<double> m_pulseFallTime_ns{this, "PulseFallTime_ns", 7.0,
          "Fall time of the single electron pulse (ns)"};

  //Analog waveform window
  Gaudi::Property<double> m_waveformTimeWindow_ns{this, "WaveformTimeWindow_ns", 300.0,
          "Time window of analog waveform (ns)"};

  //Analog waveform bin size
  Gaudi::Property<double> m_waveformBinSize_ns{this, "WaveformBinSize_ns", 1.0,
          "Time bin size for the analog waveform (ns)"};

  //Function used to produce Single Pulse
  double singleElectronPulse(double t_ns, double t0_ns, double q) const;

  // --------------- Signal propagation along sense wire  ----------------------//
  
  //Signal speed along the wire
  Gaudi::Property<double> m_signalSpeed_mm_per_ns{this, "signalSpeed_mm_per_ns", 200.0,
	  "Propagation speed along sense wire [mm/ns]. 200 mm/ns ~ 5 ns/m"};

  //Attenuation in signal
  Gaudi::Property<double> m_wireAttenuationLength_mm{this, "WireAttenuationLength_mm", 1.0e6,
	  "Attenuation length lambda [mm] for exp(-d/lambda). Very large disables attenuation"};
  
  //Enable or disable the attenuation
  Gaudi::Property<bool> m_enableWireAttenuation{this, "EnableWireAttenuation", true,
	  "Enable signal attenuation along the sense wire (charge-division effect)"};

  // ---------- impedence mismatch + electronics transfer fucntion: -----------//

  //Front end gain of the electronics
  Gaudi::Property<double> m_frontEndGain{this, "FrontEndGain", 1.0,
	  "Overall analog gain applied before ADC (dimensionless)"};

  //Apply mismatch
  Gaudi::Property<bool> m_enableImpedanceMismatch{this, "EnableImpedanceMismatch", true,
	  "Enable impedance mismatch scaling (reflection/transmission)"};

  //Set the value of mismatch resistance
  Gaudi::Property<double> m_matchingRes_Ohm{this, "MatchingRes_Ohm", 330.0,
	  "Matching / termination resistance at preamp input [Ohm]"};

  //Set the value of the Tube Impedance
  Gaudi::Property<double> m_tubeImpedance_Ohm{this, "TubeImpedance_Ohm", 50.0,
	  "Transmission line / tube characteristic impedance [Ohm]"};

  //Apply Electronic transfer function
  Gaudi::Property<bool> m_enableElectronicsTF{this, "EnableElectronicsTF", true,
	  "Enable electronics transfer function (shaper impulse response)"};

  //Set the value of Time rise
  Gaudi::Property<double> m_elecTauRise_ns{this, "ElectronicsTauRise_ns", 3.0,
	  "Electronics shaping rise time constant (ns)"};

  //Set the value of Time fall1
  Gaudi::Property<double> m_elecTauFall1_ns{this, "ElectronicsTauFall1_ns", 9.0,
          "Electronics shaping fall time constant #1 (ns)"};

  //Set the value of Time fall2
  Gaudi::Property<double> m_elecTauFall2_ns{this, "ElectronicsTauFall2_ns", 25.0,
          "Electronics shaping fall time constant #2 (ns)"};

  //Set the fraction value of the Time fall1 and fall2
  Gaudi::Property<double> m_elecMixFraction{this, "ElectronicsMixFraction", 0.5,
	  "Mixing fraction between fall1 and fall2 (0..1)"};

  //
  Gaudi::Property<double> m_elecKernelLength_ns{this, "ElectronicsKernelLength_ns", 200.0,
	  "Length of the impulse response kernel used for discrete convolution (ns)"};

  //--------- FFT/IFFT colored noise (from par.root/fft_noise) --------------//
  
  Gaudi::Property<bool> m_enableFFTNoise{this, "EnableFFTNoise", true,
	  "Enable adding colored noise using FFT magnitude from par.root"};
  
  Gaudi::Property<std::string> m_noiseFileName{this, "NoiseFileName", "par.root",
	  "ROOT file that contains directory fft_noise with fft_freq/fft_mag/fft_amp"};
  
  Gaudi::Property<std::string> m_noiseDirName{this, "NoiseDirName", "fft_noise",
	  "Directory name inside ROOT file containing fft_freq/fft_mag/fft_amp"};
  
  Gaudi::Property<double> m_noiseScale{this, "NoiseScale", 1.0e-3,
	  "Overall scale factor for FFT noise (multiplies the generated time-noise)"};
  
  Gaudi::Property<bool> m_noiseRemoveDC{this, "NoiseRemoveDC", true,
	  "Remove DC offset from generated noise (recommended)"};
  // Cached FFT-noise template (loaded once in initialize)
  mutable bool m_fftNoiseReady{false};
  mutable std::vector<double> m_fftMag;
  mutable double m_fftAmpNorm{1.0};     // fft_amp (normalization from file)
  mutable double m_fftMaxFreq{0.0};     // max frequency from fft_freq (same units as 1/ns)
  mutable int    m_fftSize{0};          // size of m_fftMag

  // --- Digitization of the analog waveform:
  Gaudi::Property<bool>   m_doWaveformDigitization{this, "DoWaveformDigitization", true,
  "Digitize analog waveform into ADC time bins"};
  
  Gaudi::Property<int> m_adcPolarity{this, "ADCPolarity", +1,
  "+1: positive pulse, -1: negative pulse (electronics polarity)"};
  
  Gaudi::Property<int>    m_adcBits{this, "ADCBits",12, 
	  "ADC resolution bits"};
  
  Gaudi::Property<double> m_adcLSB_mV{this, "ADCLSB_mV", 0.5, 
	  "ADC step (LSB) in mV"};
  
  Gaudi::Property<double> m_adcBaseline_mV{this, "ADCBaseline_mV", 50.0, 
	  "Baseline in mV"};


  //-----------------------------------End of Changes--------------------------------------//
  
  //        debug information

  /// Flag to create output file with debug histgrams
  Gaudi::Property<bool> m_create_debug_histos{this, "create_debug_histograms", false,
                                              "Create output file with histograms for debugging"};

  /// name for the file that will contain the histograms for debugging
  Gaudi::Property<std::string> m_out_debug_filename{this, "out_debug_filename", "dch_digi_alg_debug.root",
                                                    "name for the file that will contain the histograms for debugging"};
  /// histogram to store distance from sim hit position to the sense wire
  TH1D* hDpw;

  /// histogram to store distance from digi-hit to the wire. Should be zero because digi-hit central position lies on
  /// the wire. This histogram is a consistency check, because the function used to calculate the distance to the wire
  /// is different from the function used to calculate the digi-hit central position from a sim-hit position
  TH1D* hDww;

  /// histogram to store smearing along the wire
  TH1D* hSz;

  /// histogram to store smearing perpendicular the wire
  TH1D* hSxy;

//------------------------Changes made by Muhammad Saiel----------------------//
  TGraph* gSenseWireXY{nullptr};
  TH1F* hPathLength;
  TH2F* hPLvsnC;
  TH1F* hX;
  TH1F* hY;
  TH1F* hZ;
  TH3F* hXYZ;
  TH1F* heDep;
  TH1F* hR;
  TH2F* hRvsZ;
  TH1F* hTotPathCell;
  TH1F* hBetaGammaCell;

  TH1F* hNcl_perStep;
  TH1F* hNe_perStep;
  TH1F* hClSpacing_mm;
  
  TH1F* hTd;
  TProfile* hXT;
  TGraph2D* m_grXT2D = nullptr;
  TH1F* hV;
  TProfile* pVvsE;
  TH2F* hVvsR;

  TH1F* hAvalancheQ;

  mutable TH1F* hSignalPulse;
  TH1F* hWaveform;
  TH1F* hWaveformL;
  TH1F* hWaveformR;
  TH1F* hWaveformElecL{nullptr};
  TH1F* hWaveformElecR{nullptr};
  TH1F* hWaveformDigiL;
  TH1F* hWaveformDigiR;
  TH1F* hADCL;
  TH1F* hADCR;

  mutable bool m_waveformFilled{false};

  mutable TTree* m_clusterTree{nullptr};
  mutable float m_cluster_Ncl_step{0.f};

//--------------------------------------End of Changes-------------------------//

  /// Create ROOT file for debug histograms
  /// Does not change ROOT directory
  void Create_outputROOTfile_for_debugHistograms();
};

DECLARE_COMPONENT(DCHdigi_v01);

#endif
