#include "DCHXT2DLUT.h"

#include "TFile.h"
#include "TGraph2D.h"
#include "TRandom3.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

DCHXT2DLUT::DCHXT2DLUT(const std::string& rootFile,
                       const std::string& meanName,
                       const std::string& sigmaName) {
  // Open ROOT file and keep ownership in a unique_ptr
  m_file.reset(TFile::Open(rootFile.c_str(), "READ"));
  if (!m_file || m_file->IsZombie()) {
    throw std::runtime_error("[DCHXT2DLUT] ERROR: cannot open file " + rootFile);
  }

  // Retrieve the 2D x-t mean and sigma surfaces
  m_file->GetObject(meanName.c_str(), m_gMean);
  m_file->GetObject(sigmaName.c_str(), m_gSigma);

  if (!m_gMean || !m_gSigma) {
    throw std::runtime_error("[DCHXT2DLUT] ERROR: missing TGraph2D '" + meanName +
                             "' or '" + sigmaName + "' in " + rootFile);
  }

  std::cout << "[DCHXT2DLUT] Loaded " << meanName << " and " << sigmaName
            << " from " << rootFile << "\n";
}

/// Helper function: checks if the LUT query result looks invalid.
bool DCHXT2DLUT::isInvalid(double mu, double sig, double x, double y) {
  if (!std::isfinite(mu) || !std::isfinite(sig)) return true;
  if (mu < 0.0 || sig < 0.0) return true;

  // Typical TGraph2D behavior outside the interpolation region:
  // values can collapse to ~0, so reject that away from the wire center.
  const double r = std::hypot(x, y);
  if (r > 1e-4 && std::abs(mu) < 1e-12 && std::abs(sig) < 1e-12) return true;

  return false;
}

/// Query mean and sigma drift time at a given (x,y) point in cm.
bool DCHXT2DLUT::meanSigma(double x_cm, double y_cm, double& mu_ns, double& sig_ns) const {
  if (!m_gMean || !m_gSigma) {
    mu_ns = 0.0;
    sig_ns = 0.0;
    return false;
  }

  // Interpolate the two lookup surfaces at the requested transverse position
  mu_ns = m_gMean->Interpolate(x_cm, y_cm);
  sig_ns = m_gSigma->Interpolate(x_cm, y_cm);

  // Reject non-physical or invalid interpolation results
  if (isInvalid(mu_ns, sig_ns, x_cm, y_cm)) {
    mu_ns = 0.0;
    sig_ns = 0.0;
    return false;
  }

  return true;
}

/// Sample drift time in ns using the digitizer RNG.
/// 1) Read mean and sigma from the LUT
/// 2) If sigma > 0, apply Gaussian smearing
/// 3) If sigma == 0, use the mean directly
/// 4) Enforce non-negative drift time
double DCHXT2DLUT::sampleTimeNs(double x_cm, double y_cm, TRandom3& rng) const {
  double mu = 0.0;
  double sig = 0.0;

  const bool ok = meanSigma(x_cm, y_cm, mu, sig);
  if (!ok) return 0.0;

  double t_ns = (sig > 0.0) ? rng.Gaus(mu, sig) : mu;
  if (t_ns < 0.0) t_ns = 0.0;

  return t_ns;
}

