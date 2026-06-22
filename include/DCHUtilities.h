// @author: Muhammad saiel
 
#ifndef DCH_PHYSICS_TOOLS_H
#define DCH_PHYSICS_TOOLS_H
#include <vector>
//#include <random>
#include <numeric>
#include <cmath>
#include <edm4hep/MCParticle.h>
#include <array>
#include "TRandom3.h"
#include <algorithm>

// Function compute betaGamma:
inline double compute_beta_gamma(const edm4hep::MCParticle& mc)
{
  const auto& p = mc.getMomentum();        // GeV
  const double p_mag = std::sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
  const double m     = mc.getMass();       // GeV
  if (m <= 0.0) return 0.0;
  return p_mag / m;
}

namespace
{

inline int ZT_poisson(double mu, TRandom3& rng)
{
  if (!(mu > 0) || !std::isfinite(mu)) return 0;
  
  int k = 0;
  do {
    k = rng.Poisson(mu); 
  }  while (k == 0);
  return k;
}

// Generate exponentially-distributed cluster positions along a step of length l_mm.
// Output: vector of positions in mm, measured from the start of the step.
inline std::vector<double>
generate_cluster_positions_mm(double l_mm, int nClusters, TRandom3& rng)
{
  std::vector<double> pos_mm;
  if (l_mm <= 0.0 || nClusters <= 0) return pos_mm;

  pos_mm.reserve(nClusters);

  // Conditional on having exactly nClusters in the step,
  // their positions are the sorted order statistics of Uniform(0, l_mm).
  for (int i = 0; i < nClusters; ++i) {
    pos_mm.push_back(rng.Uniform(0.0, l_mm));
  }

  std::sort(pos_mm.begin(), pos_mm.end());
  return pos_mm;
}

// --- Experimental cluster size probabilities w(n) for He-iSobutane (90:10) ---
// Based on Fischle et al., NIM A301 (1991)
inline constexpr std::array<double, 35> make_w_cluster_He_iC4H10()
{
  std::array<double, 35> w{};

  // Tabulated values for n = 1 ... 19 (given in %)
  w[0]  = 76.60;
  w[1]  = 12.50;
  w[2]  =  4.60;
  w[3]  =  2.00;
  w[4]  =  1.20;
  w[5]  =  0.75;
  w[6]  =  0.50;
  w[7]  =  0.36;
  w[8]  =  0.25;
  w[9]  =  0.19;
  w[10] =  0.14;
  w[11] =  0.10;
  w[12] =  0.08;
  w[13] =  0.06;
  w[14] =  0.048;
  w[15] =  0.043;
  w[16] =  0.038;
  w[17] =  0.034;
  w[18] =  0.030;

  // Tail formula for n >= 20 : 10.9 / n^2  (still in % here)
  for (int k = 20; k <= 35; ++k) {
    w[k - 1] = 10.9 / (double(k) * double(k));
  }

  // Convert % -> fraction
  for (auto& x : w) {
    x *= 0.01;
  }

  return w;
}

inline constexpr auto w_cluster_He_iC4H10 = make_w_cluster_He_iC4H10();
//----------------------------------------------------------------------------------//
inline double mean_cluster_size_He()
{
  double Ne_mean = 0.0;
  for (size_t i = 0; i < w_cluster_He_iC4H10.size(); ++i) {
    Ne_mean += (double(i) + 1.0) * w_cluster_He_iC4H10[i];
  }
  return Ne_mean;
}

// --- Draw a cluster size n_e according to experimental probabilities ---
inline int sample_cluster_size(TRandom3& rng)
{
  double total = 0.0;
  for (double w : w_cluster_He_iC4H10) {
    total += w;
  }

  if (!(total > 0.0)) {
    return 1;
  }

  const double u = rng.Rndm() * total;

  double cdf = 0.0;
  for (size_t i = 0; i < w_cluster_He_iC4H10.size(); ++i) {
    cdf += w_cluster_He_iC4H10[i];
    if (u <= cdf) {
      return static_cast<int>(i) + 1;
    }
  }

  return static_cast<int>(w_cluster_He_iC4H10.size());
}

} //end namespace

#endif
