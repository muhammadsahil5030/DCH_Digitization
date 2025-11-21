// @ auther: Muhammad saiel
//
#ifndef DCH_PHYSICS_TOOLS_H
#define DCH_PHYSICS_TOOLS_H

#include <vector>
#include <random>
#include <numeric>
#include <cmath>

//include EDM4hep
#include <edm4hep/MCParticle.h>

struct CellAcc
{
  double mu_sum = 0.0;
  double length_mm = 0.0;
  double bgL_sum = 0.0;
};

// Function compute betaGamma:
inline double compute_beta_gamma(const edm4hep::MCParticle& mc)
{
  const auto& p = mc.getMomentum();        // GeV
  const double p_mag = std::sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
  const double m     = mc.getMass();       // GeV
  if (m <= 0.0) return 0.0;
  return p_mag / m;
}

// This function is used to calculate the dNdx using BB.h file:
// when we have a beam-test TGraph gNcldx(bg)->(clusters/cm),
double DCHdigi_v01::get_dNcldx_per_cm(double betagamma) const
{
  // 1) reconstruct momentum from β
  const double m_GeV = m_MassForBB_GeV.value();
  const double p_GeV = betagamma * m_GeV;

  // 2) call your header: bethe_bloch(xp, par)
  double xp[1]  = { p_GeV };                           // GeV (header multiplies by 1e3 internally)
  double par[2] = { m_GeV * 1.0e3,                     // mass in MeV
                    m_MeanExcEnergy_eV.value() * 1e-6  // I in MeV
                  };
  const double dEdx_MeVcm2_per_g = BB::bethe_bloch(xp, par);
  const double dEdx_MeV_per_cm = dEdx_MeVcm2_per_g * m_GasDensity_g_cm3.value();
  const double lambda_per_cm = (dEdx_MeV_per_cm * 1.0e6) / m_W_eff_eV.value();

  return (lambda_per_cm > 0.0) ? lambda_per_cm : 0.0;
}

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
  std::exponential_distribution<double> expo(lambda_per_cm);

  double s_cm = 0.0;
  while (true) {
    s_cm += expo(rng);           // draw next spacing
    if (s_cm >= l_cm) break;     // stop if we exceed the step length
    pos_mm.push_back(10.0 * s_cm); // store position in mm
  }
  return pos_mm;

}

// --- Experimental cluster size probabilities w(n) for He-iSobutane (90:10) ---
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

#endif
