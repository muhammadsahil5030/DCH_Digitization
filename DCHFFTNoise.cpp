#include "DCHFFTNoise.h"

#include <algorithm>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TVirtualFFT.h"

namespace dchfft {

static inline int float_cmp(double a, double b, double eps=1e-9) {
  if (std::fabs(a-b) < eps) return 0;
  return (a > b) ? 1 : -1;
}

bool loadFFTNoiseTemplate(
    const std::string& rootFile,
    const std::string& treeName,
    std::vector<double>& magTemplate,
    double& ampNorm,
    double& maxFreq,
    int& size,
    std::string* errMsg
) {
  magTemplate.clear();
  ampNorm = 1.0;
  maxFreq = 0.0;
  size    = 0;

  TFile f(rootFile.c_str(), "READ");
  if (f.IsZombie()) {
    if (errMsg) *errMsg = "Cannot open ROOT file: " + rootFile;
    return false;
  }

  TTree* t = dynamic_cast<TTree*>(f.Get(treeName.c_str()));
  if (!t) {
    if (errMsg) *errMsg = "Cannot find TTree '" + treeName + "' in " + rootFile;
    return false;
  }

  std::vector<double>* fft_freq = nullptr;
  std::vector<double>* fft_mag  = nullptr;
  double fft_amp = 0.0;

  t->SetBranchAddress("fft_freq", &fft_freq);
  t->SetBranchAddress("fft_mag",  &fft_mag);
  t->SetBranchAddress("fft_amp",  &fft_amp);

  if (t->GetEntries() <= 0) {
    if (errMsg) *errMsg = "TTree '" + treeName + "' has 0 entries.";
    return false;
  }

  t->GetEntry(0);

  if (!fft_freq || !fft_mag || fft_freq->empty() || fft_mag->empty() || fft_amp == 0.0) {
    if (errMsg) *errMsg = "Bad branches: need fft_freq, fft_mag (non-empty) and fft_amp != 0";
    return false;
  }

  magTemplate = *fft_mag;
  size        = (int)magTemplate.size();
  maxFreq     = fft_freq->back();
  ampNorm     = fft_amp;

  return true;
}

std::vector<double> makeFFTNoiseSamples(
    int nSamples,
    double dt_ns,
    const std::vector<double>& magTemplate,
    double maxFreqTemplate,
    TRandom3& rnd,
    bool removeDC
) {
  std::vector<double> out;
  out.reserve(nSamples);

  // waveform sampling frequency and Nyquist
  const double f_wf = 1.0 / dt_ns;      // 1/ns
  const double f_ny = 0.5 * f_wf;       // 1/ns

  // local working copy (do not modify cached template)
  int fftSize = (int)magTemplate.size();
  double fftMaxFreq = maxFreqTemplate;
  int downFactor = 1;

  std::vector<double> mag = magTemplate;

  // Match template freq range to waveform Nyquist (same logic as NoiseGenerator.C)
  const int cmp = float_cmp(fftMaxFreq, f_ny);
  if (cmp < 0) {
    // template covers less freq -> pad with zeros
    fftSize = (int)(f_ny / fftMaxFreq * fftSize);
    fftMaxFreq = f_ny;
    mag.resize(fftSize, 0.0);
  } else if (cmp > 0) {
    // template covers higher freq -> upsample then downsample
    int scale_factor = (int)(fftMaxFreq / f_wf * 2.0) + 1;
    fftSize = (int)(f_ny * scale_factor / fftMaxFreq * fftSize);
    fftMaxFreq = f_ny * scale_factor;
    downFactor = scale_factor;
    mag.resize(fftSize, 0.0);
  }

  // time window covered by one IFFT segment
  const double dt_seg = 1.0 / (2.0 * fftMaxFreq);
  const double segWindow_ns = fftSize * dt_seg;

  const double wfWindow_ns = nSamples * dt_ns;
  const int nSegments = (segWindow_ns < wfWindow_ns) ? (int)(wfWindow_ns / segWindow_ns) + 1 : 1;

  std::vector<double> tmp;
  tmp.reserve((size_t)nSegments * (size_t)fftSize);

  for (int iseg = 0; iseg < nSegments; ++iseg) {
    TVirtualFFT* ifft = TVirtualFFT::FFT(1, &fftSize, "C2R M K");

    std::vector<double> re(fftSize, 0.0);
    std::vector<double> im(fftSize, 0.0);

    for (int i = 0; i < fftSize; ++i) {
      const double ph = rnd.Rndm() * 2.0 * TMath::Pi();
      re[i] = mag[i] * std::cos(ph);
      im[i] = mag[i] * std::sin(ph);
    }

    ifft->SetPointsComplex(&re[0], &im[0]);
    ifft->Transform();
    ifft->GetPointsComplex(&re[0], &im[0]); // re now holds time samples

    tmp.insert(tmp.end(), re.begin(), re.end());
    delete ifft;
  }

  // Downsample if needed (average blocks)
  if (downFactor > 1) {
    std::vector<double> ds;
    ds.reserve(tmp.size() / downFactor + 1);
    for (size_t i = 0; i < tmp.size(); i += (size_t)downFactor) {
      double avg = 0.0;
      int cnt = 0;
      for (int j = 0; j < downFactor && (i + (size_t)j) < tmp.size(); ++j) {
        avg += tmp[i + (size_t)j];
        ++cnt;
      }
      avg /= std::max(1, cnt);
      ds.push_back(avg);
    }
    tmp.swap(ds);
  }

  out.assign(tmp.begin(), tmp.begin() + std::min((int)tmp.size(), nSamples));
  if ((int)out.size() < nSamples) out.resize(nSamples, 0.0);

  if (removeDC && !out.empty()) {
    double mean = 0.0;
    for (double v : out) mean += v;
    mean /= (double)out.size();
    for (double& v : out) v -= mean;
  }

  return out;
}

} // namespace dchfft

