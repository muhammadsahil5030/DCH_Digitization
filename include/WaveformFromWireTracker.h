/** ================= WaveformFromWireTracker =================
 * Gaudi Algorithm for wire-tracker waveform production.
 *
 * @authors: Muhammad Saiel, Giovanni Francesco Tassielli, Nicola De Filippis
 * @date   2026-04
 *
 * <h4>Purpose</h4>
 * This algorithm is the waveform-producing stage of the digitization chain.
 *
 * It consumes gas-ionization clusters produced upstream by
 * DistributeClustersInGas and builds the detector response seen at the two
 * readout ends of each sense wire.
 *
 * Each input object represents one generated gas-ionization cluster. For each
 * cluster, the algorithm computes the drift time, samples the electron
 * multiplication, builds the corresponding signal contribution, and accumulates
 * all contributions belonging to the same cell.
 *
 * The accumulated cell signal is propagated to the left and right ends of the
 * sense wire, shaped by the front-end electronics response, combined with
 * FFT-based colored electronic noise, and digitized into ADC samples.
 *
 * The final output is a pair of edm4hep::RawTimeSeriesCollection objects:
 * one collection for the left-end waveform and one collection for the
 * right-end waveform.
 *
 * This code uses DD4hep natural length units (cm), while EDM4hep positions and
 * distances are stored in mm. Unit conversions must therefore be handled
 * explicitly.
 *
 * <h4>Input collections and prerequisites</h4>
 * Processor requires:
 *   - edm4hep::SimTrackerHitCollection produced upstream by
 *     DistributeClustersInGas, where each SimTrackerHit represents one
 *     generated gas-ionization cluster.
 *   - edm4hep::EventHeaderCollection used to build reproducible event-level
 *     random seeds.
 *
 * <h4>Scope of this stage</h4>
 * This algorithm is responsible for:
 *   - reading gas-ionization clusters;
 *   - computing drift times from the x-t lookup table;
 *   - sampling avalanche charge fluctuations with a Polya model;
 *   - accumulating all drifted electrons belonging to the same readout cell;
 *   - building analog waveforms from individual electron pulses;
 *   - propagating the signal to the left and right wire ends;
 *   - applying front-end/electronics response;
 *   - adding FFT-based colored electronic noise;
 *   - digitizing the waveform into ADC samples;
 *   - producing left/right edm4hep::RawTimeSeriesCollection outputs.
 *
 * <h4>Output</h4>
 * Processor produces:
 *   - edm4hep::RawTimeSeriesCollection for the digitized left-end waveform;
 *   - edm4hep::RawTimeSeriesCollection for the digitized right-end waveform.
 *
 * <h4>Main input/output collection properties</h4>
 * @param InputClusterCollection Input gas-cluster collection,
 * type edm4hep::SimTrackerHitCollection.
 *
 * @param DigiWaveformLCollection Output digitized waveform collection for
 * the left wire end, type edm4hep::RawTimeSeriesCollection.
 *
 * @param DigiWaveformRCollection Output digitized waveform collection for
 * the right wire end, type edm4hep::RawTimeSeriesCollection.
 *
 * <h4>Main configurable properties</h4>
 * @param DetectorName Detector/subdetector name in DD4hep.
 *
 * @param XTFileName ROOT file containing the x-t lookup table used by
 * DCHXT2DLUT.
 *
 * @param GasGain Mean gas gain used in the Polya avalanche model.
 *
 * @param PolyaTheta Shape parameter of the Polya avalanche model.
 *
 * @param PulseAmplitudeScale Global scale converting avalanche charge to
 * voltage amplitude.
 *
 * @param WaveformStartTime_ns Start time of the waveform time window in ns.
 *
 * @param WaveformTimeWindow_ns Total waveform time window in ns.
 *
 * @param WaveformBinSize_ns Analog waveform bin size in ns.
 *
 * @param EnableFFTNoise Enable FFT-based colored noise injection.
 *
 * @param NoiseFileName ROOT file containing the FFT noise template.
 *
 * @param NoiseDirName Directory/key name for FFT noise objects in the ROOT file.
 *
 * @param DoWaveformDigitization Enable conversion of analog waveform samples
 * into ADC samples.
 *
 * @param ADCSamplePeriod_ns ADC sampling period in ns.
 *
 * @param ADCPolarity Electronics polarity convention.
 *
 * @param ADCBits ADC resolution in bits.
 *
 * @param ADCLSB_mV ADC least-significant-bit size in mV.
 *
 * @param ADCBaseline_mV ADC baseline offset in mV.
 *
 * @param create_debug_histograms Optional flag to enable QA/debug histograms
 * registered to THistSvc. The output ROOT file name is configured from the
 * Gaudi steering file.
 *
 * @param GeoSvcName Geometry service name.
 *
 * @param uidSvcName UniqueIDGenSvc instance name used for reproducible seeds.
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
#include "edm4hep/RawTimeSeriesCollection.h"
#include "edm4hep/utils/vector_utils.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DDSegmentation/BitFieldCoder.h"


// DD4hep detector extension
#include "DDRec/DCH_info.h"

// Drift time lookup table
#include "DCHXT2DLUT.h"

// ROOT headers
#include "TH1F.h"
#include "TRandom3.h"
#include "TVector3.h"

// STL
#include <memory>
#include <string>
#include <vector>

class TF1;

struct WaveformFromWireTracker final
    : k4FWCore::MultiTransformer<
          std::tuple<
	  // Digitized waveform collections at the two wire ends
	  edm4hep::RawTimeSeriesCollection,
	  edm4hep::RawTimeSeriesCollection>(
			  const edm4hep::SimTrackerHitCollection&,
			  const edm4hep::EventHeaderCollection&)> {

  WaveformFromWireTracker(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  StatusCode finalize() override;

  std::tuple<
	  //Digi waveform collections
	  edm4hep::RawTimeSeriesCollection,
	  edm4hep::RawTimeSeriesCollection>
  operator()(const edm4hep::SimTrackerHitCollection&,
		  const edm4hep::EventHeaderCollection&) const override;

private:
  /// conversion factor mm to cm, static to the class to avoid clash with DD4hep
  static constexpr double MM_TO_CM = 0.1;

  //------------------------------------------------------------------
  //          machinery for geometry

  /// Geometry service name
  Gaudi::Property<std::string> m_geoSvcName{this, "GeoSvcName", "GeoSvc", "The name of the GeoSvc instance"};
  Gaudi::Property<std::string> m_uidSvcName{this, "uidSvcName", "uidSvc", "The name of the UniqueIDGenSvc instance"};

  /// Detector name
  Gaudi::Property<std::string> m_detectorName{this, "DetectorName", "DCH_v2", "Name of the detector"};

  /// Pointer to the geometry service
  SmartIF<IGeoSvc> m_geoSvc;

  /// Decoder for the cellID
  const dd4hep::DDSegmentation::BitFieldCoder* m_decoder{nullptr};

  /// Pointer to drift chamber data extension
  dd4hep::rec::DCH_info* dch_data = {nullptr};

  double m_halfChamberLength_mm = 0.0;	// gas half-length along z

  /// create seed using the uid
  SmartIF<IUniqueIDGenSvc> m_uidSvc;

  /// Create random engine, initialized with seed out of Event Header
  TRandom3 CreateRandomEngine(const edm4hep::EventHeaderCollection& headers) const;
  //------------------------------------------------------------------

  /// Print algorithm configuration
  void PrintConfiguration(std::ostream& io) const;

  int CalculateLayerFromCellID(dd4hep::DDSegmentation::CellID id) const {
    return dch_data->CalculateILayerFromCellIDFields(m_decoder->get(id, "layer"), m_decoder->get(id, "superlayer"));
  }

  int CalculateNphiFromCellID(dd4hep::DDSegmentation::CellID id) const { return m_decoder->get(id, "nphi"); }

  TVector3 Convert_EDM4hepVector_to_TVector3(const edm4hep::Vector3d& v, double scale) const {
    return TVector3(v.x * scale, v.y * scale, v.z * scale);
  }

  //// ROOT file containing the x-t relation
  Gaudi::Property<std::string> m_xtFileName{this, "XTFileName", "par.root",
	  "ROOT file with x-t relation used to convert radius to drift time" };

  //pointer to the LUT
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

  //------------------------- pulse shaping block -------------------------//

  //Scale parameter
  Gaudi::Property<double> m_pulseAmplitudeScale{this, "PulseAmplitudeScale", 5.e-6,
          "Global Scale factor converting avalache charge to voltage"};
  
  //Pulse rise time
  Gaudi::Property<double> m_pulseRiseTime_ns{this, "PulseRiseTime_ns", 1.0,
          "Rise time of the single electron pulse (ns)"};

  //Pulse decay time
  Gaudi::Property<double> m_pulseFallTime_ns{this, "PulseFallTime_ns", 7.0,
          "Fall time of the single electron pulse (ns)"};

  //Pulse bin size
  Gaudi::Property<double> m_pulseBinSize_ns{this, "PulseBinSize_ns", 1.0,
          "Time bin size for the analog pulse (ns)"};

    // Time window used only for the debug single-electron pulse histogram
  Gaudi::Property<double> m_signalPulseWindow_ns{this, "SignalPulseWindow_ns", 500.0,
          "Displayed time window of the debug single-electron pulse histogram (ns)"};

  // Analog waveform start time
  Gaudi::Property<double> m_waveformStartTime_ns{this, "WaveformStartTime_ns", 0.0,
	  "Start time of the waveform window in ns" };

  //Analog waveform window
  Gaudi::Property<double> m_waveformTimeWindow_ns{this, "WaveformTimeWindow_ns", 3000.0,
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
  Gaudi::Property<double> m_wireAttenuationLength_mm{this, "WireAttenuationLength_mm", 1.0e4,
	  "Attenuation length lambda [mm] for exp(-d/lambda). Very large disables attenuation"};

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

  // Length (ns) of the electronics impulse response used in waveform convolution.
  Gaudi::Property<double> m_elecKernelLength_ns{this, "ElectronicsKernelLength_ns", 200.0,
	  "Length of the impulse response kernel used for discrete convolution (ns)"};

  //--------- FFT/IFFT colored noise (from par.root/fft_noise) --------------//
  
  // Enable addition of frequency-shaped (colored) electronic noise using FFT/IFFT.
  Gaudi::Property<bool> m_enableFFTNoise{this, "EnableFFTNoise", true,
	  "Enable adding colored noise using FFT magnitude from par.root"};
  
  // ROOT file containing the FFT-based noise frequency template.
  Gaudi::Property<std::string> m_noiseFileName{this, "NoiseFileName", "par.root",
	  "ROOT file that contains directory fft_noise with fft_freq/fft_mag/fft_amp"};
  
  // Directory name inside the ROOT file holding the FFT noise histograms.
  Gaudi::Property<std::string> m_noiseDirName{this, "NoiseDirName", "fft_noise",
	  "Directory name inside ROOT file containing fft_freq/fft_mag/fft_amp"};
  
  // Global scaling factor applied to the generated time-domain noise amplitude.
  Gaudi::Property<double> m_noiseScale{this, "NoiseScale", 1.0e-3,
	  "Overall scale factor for FFT noise (multiplies the generated time-noise)"};
  
  // Remove DC (baseline) component from the generated noise waveform.
  Gaudi::Property<bool> m_noiseRemoveDC{this, "NoiseRemoveDC", true,
	  "Remove DC offset from generated noise (recommended)"};
  // Cached FFT-noise template (loaded once in initialize)
  bool m_fftNoiseReady{false};
  std::vector<double> m_fftMag;
  double m_fftAmpNorm{1.0};     // fft_amp (normalization from file)
  double m_fftMaxFreq{0.0};     // max frequency from fft_freq (same units as 1/ns)
  int    m_fftSize{0};          // size of m_fftMag

  // Enable conversion of the analog waveform into a digitized ADC waveform.
  Gaudi::Property<bool>   m_doWaveformDigitization{this, "DoWaveformDigitization", true,
  "Digitize analog waveform into ADC time bins"};

  // ADC sampling period for digitized waveform
  Gaudi::Property<double> m_adcSamplePeriod_ns{this, "ADCSamplePeriod_ns", 0.5,
	  "ADC sampling period (ns) for digitized waveform"};
  
  // Set electronics polarity: +1 keeps pulse sign, -1 inverts it.
  Gaudi::Property<int> m_adcPolarity{this, "ADCPolarity", +1,
  "+1: positive pulse, -1: negative pulse (electronics polarity)"};
  
  // Number of bits defining the ADC amplitude resolution.
  Gaudi::Property<int>    m_adcBits{this, "ADCBits",12, 
	  "ADC resolution bits"};
  
  // Voltage step size (mV) corresponding to one ADC count
  Gaudi::Property<double> m_adcLSB_mV{this, "ADCLSB_mV", 0.5, 
	  "ADC step (LSB) in mV"};
  
  // Baseline voltage (mV) added to the digitized waveform
  Gaudi::Property<double> m_adcBaseline_mV{this, "ADCBaseline_mV", 50.0, 
	  "Baseline in mV"};

  // Enable diagnostic histograms registered to THistSvc
  Gaudi::Property<bool> m_create_debug_histos{this, "create_debug_histograms", false,
                                              "Enable diagnostic histograms via THistSvc"};

  /// histogram to store distance from sim hit position to the sense wire
  SmartIF<ITHistSvc> m_histSvc;
  TH1F* hAvalancheQ{nullptr};
  mutable TH1F* hSignalPulse{nullptr};
  TH1F* hWaveform{nullptr};
  TH1F* hWaveformL{nullptr};
  TH1F* hWaveformR{nullptr};
  TH1F* hWaveformElecL{nullptr};
  TH1F* hWaveformElecR{nullptr};
  TH1F* hNoise{nullptr};
  TH1F* hWaveformDigiL{nullptr};
  TH1F* hWaveformDigiR{nullptr};

  // Flag used only for debug histogram.
  // Once one representative waveform has been stored, the detailed waveform
  // histograms are not overwritten by later cells or later events, 
  // to avoid overfilling and to preserve a clean example for debugging.
  mutable bool m_waveformFilled{false};

};

DECLARE_COMPONENT(WaveformFromWireTracker);


