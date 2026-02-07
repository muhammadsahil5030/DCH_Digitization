#pragma once
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

class TRandom3;  // forward declare

class DCHXT2DLUT {
public:
  bool load(const std::string& rootFile,
            const std::string& meanName  = "xt_mean",
            const std::string& sigmaName = "xt_error",
            int nx = 240, int ny = 240);

  bool isLoaded() const { return m_loaded; }

  bool meanSigma(double x_cm, double y_cm, double& mu_ns, double& sig_ns) const;

  // use same RNG as DCHdigi
  double sampleTimeNs(double x_cm, double y_cm, TRandom3& rng) const;

  double rmaxCm() const { return m_rmax_cm; }

private:
  inline int idx(int ix, int iy) const { return iy * m_nx + ix; }

  bool bilinearAt(double x_cm, double y_cm, double& mu_ns, double& sig_ns) const;
  bool radialMeanSigma(double r_cm, double& mu_ns, double& sig_ns) const;
  void buildRadialTable(int nr = 400, int nphi = 180);

  bool m_loaded = false;
  int m_nx = 0, m_ny = 0;
  double m_xmin_cm = 0.0, m_xmax_cm = 0.0;
  double m_ymin_cm = 0.0, m_ymax_cm = 0.0;
  double m_dx = 0.0, m_dy = 0.0;
  double m_rmax_cm = 0.0;

  std::vector<float> m_mu;   // ns
  std::vector<float> m_sig;  // ns

  // 1D radial representation built from the 2D grid:
  std::vector<float> m_rbin_cm;
  std::vector<float> m_mu_r_ns;
  std::vector<float> m_sig_r_ns;
};

