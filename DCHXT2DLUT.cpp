#include "DCHXT2DLUT.h"

#include "TFile.h"
#include "TGraph2D.h"

#include <iostream>
#include <limits>
#include <cmath>

#include "TRandom3.h"

namespace {

// treat ROOT out-of-domain interpolation (often returns ~0) as invalid
static inline bool isInvalid(double mu, double sg, double x, double y) {
  if (!std::isfinite(mu) || !std::isfinite(sg)) return true;
  if (mu < 0.0 || sg < 0.0) return true;
  const double r = std::hypot(x, y);
  if (r > 1e-4 && std::abs(mu) < 1e-9 && std::abs(sg) < 1e-9) return true;
  return false;
}

// quadratic through 3 points (Lagrange form)
static inline double quad3(double x,
                           double x1, double x2, double x3,
                           double y1, double y2, double y3)
{
  const double L1 = (x - x2) * (x - x3) / ((x1 - x2) * (x1 - x3));
  const double L2 = (x - x1) * (x - x3) / ((x2 - x1) * (x2 - x3));
  const double L3 = (x - x1) * (x - x2) / ((x3 - x1) * (x3 - x2));
  return y1 * L1 + y2 * L2 + y3 * L3;
}

// Query the original graphs ONLY during grid build:
// - If inside convex hull => use interpolation
// - If outside => DO NOT walk inward here anymore; we let the radial table handle extrapolation later
static bool queryGraphMeanSigma_NoInwardWalk(TGraph2D& gMean,
                                            TGraph2D& gSig,
                                            double x_cm, double y_cm,
                                            double& mu_ns, double& sig_ns)
{
  mu_ns  = gMean.Interpolate(x_cm, y_cm);
  sig_ns = gSig .Interpolate(x_cm, y_cm);

  if (isInvalid(mu_ns, sig_ns, x_cm, y_cm)) {
    mu_ns = 0.0;
    sig_ns = 0.0;
    return false;
  }
  return true;
}

} // namespace


bool DCHXT2DLUT::load(const std::string& rootFile,
                      const std::string& meanName,
                      const std::string& sigmaName,
                      int nx, int ny)
{
  m_loaded = false;
  m_mu.clear();
  m_sig.clear();
  m_rbin_cm.clear();
  m_mu_r_ns.clear();
  m_sig_r_ns.clear();

  if (nx < 10 || ny < 10) return false;

  TFile f(rootFile.c_str(), "READ");
  if (f.IsZombie()) return false;

  TGraph2D* gM = nullptr;
  TGraph2D* gS = nullptr;
  f.GetObject(meanName.c_str(),  gM);
  f.GetObject(sigmaName.c_str(), gS);
  if (!gM || !gS) return false;

  // bounds from points
  const int n = gM->GetN();
  if (n <= 0) return false;

  const double* xs = gM->GetX();
  const double* ys = gM->GetY();

  double xmin = +std::numeric_limits<double>::infinity();
  double xmax = -std::numeric_limits<double>::infinity();
  double ymin = +std::numeric_limits<double>::infinity();
  double ymax = -std::numeric_limits<double>::infinity();
  double rmax = 0.0;

  for (int i = 0; i < n; ++i) {
    xmin = std::min(xmin, xs[i]);
    xmax = std::max(xmax, xs[i]);
    ymin = std::min(ymin, ys[i]);
    ymax = std::max(ymax, ys[i]);
    rmax = std::max(rmax, std::hypot(xs[i], ys[i]));
  }

  const double pad = 1e-3; // cm
  m_xmin_cm = xmin - pad;
  m_xmax_cm = xmax + pad;
  m_ymin_cm = ymin - pad;
  m_ymax_cm = ymax + pad;
  m_rmax_cm = rmax;

  m_nx = nx;
  m_ny = ny;
  m_dx = (m_xmax_cm - m_xmin_cm) / (m_nx - 1);
  m_dy = (m_ymax_cm - m_ymin_cm) / (m_ny - 1);

  m_mu.assign(m_nx * m_ny, 0.0f);
  m_sig.assign(m_nx * m_ny, 0.0f);

  // Build 2D grid:
  // We do NOT “walk inward” anymore during build.
  // Any out-of-domain cells remain 0 and will be handled by the radial table + extrapolation.
  for (int iy = 0; iy < m_ny; ++iy) {
    const double y = m_ymin_cm + iy * m_dy;
    for (int ix = 0; ix < m_nx; ++ix) {
      const double x = m_xmin_cm + ix * m_dx;

      double mu = 0.0, sg = 0.0;
      const bool ok = queryGraphMeanSigma_NoInwardWalk(*gM, *gS, x, y, mu, sg);

      if (!ok) { mu = 0.0; sg = 0.0; }
      if (!std::isfinite(mu) || mu < 0.0) mu = 0.0;
      if (!std::isfinite(sg) || sg < 0.0) sg = 0.0;

      m_mu[idx(ix, iy)]  = static_cast<float>(mu);
      m_sig[idx(ix, iy)] = static_cast<float>(sg);
    }
  }

  // Build a 1D radial table from the 2D grid (no assumptions about a single direction):
  buildRadialTable(400, 180);

  m_loaded = true;
  return true;
}


// Bilinear interpolation on the internal grid WITHOUT clamping.
// Returns false if (x,y) is outside the built grid rectangle.
bool DCHXT2DLUT::bilinearAt(double x_cm, double y_cm, double& mu_ns, double& sig_ns) const
{
  if (!m_loaded && (m_mu.empty() || m_sig.empty())) return false;

  if (x_cm < m_xmin_cm || x_cm > m_xmax_cm) return false;
  if (y_cm < m_ymin_cm || y_cm > m_ymax_cm) return false;

  const double fx = (x_cm - m_xmin_cm) / m_dx;
  const double fy = (y_cm - m_ymin_cm) / m_dy;

  int ix = static_cast<int>(std::floor(fx));
  int iy = static_cast<int>(std::floor(fy));

  if (ix < 0) ix = 0;
  if (iy < 0) iy = 0;
  if (ix >= m_nx - 1) ix = m_nx - 2;
  if (iy >= m_ny - 1) iy = m_ny - 2;

  const double tx = fx - ix;
  const double ty = fy - iy;

  const float mu00 = m_mu[idx(ix,     iy    )];
  const float mu10 = m_mu[idx(ix + 1, iy    )];
  const float mu01 = m_mu[idx(ix,     iy + 1)];
  const float mu11 = m_mu[idx(ix + 1, iy + 1)];

  const float sg00 = m_sig[idx(ix,     iy    )];
  const float sg10 = m_sig[idx(ix + 1, iy    )];
  const float sg01 = m_sig[idx(ix,     iy + 1)];
  const float sg11 = m_sig[idx(ix + 1, iy + 1)];

  const double mu0 = (1.0 - tx) * mu00 + tx * mu10;
  const double mu1 = (1.0 - tx) * mu01 + tx * mu11;
  mu_ns = (1.0 - ty) * mu0 + ty * mu1;

  const double sg0 = (1.0 - tx) * sg00 + tx * sg10;
  const double sg1 = (1.0 - tx) * sg01 + tx * sg11;
  sig_ns = (1.0 - ty) * sg0 + ty * sg1;

  if (!std::isfinite(mu_ns) || mu_ns < 0.0) mu_ns = 0.0;
  if (!std::isfinite(sig_ns) || sig_ns < 0.0) sig_ns = 0.0;

  return true;
}


// Build radial mean/sigma vs r by sampling MANY angles at each r.
// This avoids “use only y=0” type shortcuts.
void DCHXT2DLUT::buildRadialTable(int nr, int nphi)
{
  m_rbin_cm.assign(nr, 0.0f);
  m_mu_r_ns.assign(nr, 0.0f);
  m_sig_r_ns.assign(nr, 0.0f);

  if (nr < 10) return;
  if (nphi < 8) nphi = 8;

  const double rmax = m_rmax_cm;
  const double dr = (rmax > 0.0) ? (rmax / nr) : 0.0;

  double last_mu = 0.0;
  double last_sg = 0.0;

  for (int ir = 0; ir < nr; ++ir) {
    const double r = (ir + 0.5) * dr;
    m_rbin_cm[ir] = static_cast<float>(r);

    double sum_mu = 0.0;
    double sum_sg = 0.0;
    int    cnt    = 0;

    for (int ip = 0; ip < nphi; ++ip) {
      const double phi = (2.0 * M_PI) * (double(ip) / double(nphi));
      const double x = r * std::cos(phi);
      const double y = r * std::sin(phi);

      double mu = 0.0, sg = 0.0;
      if (!bilinearAt(x, y, mu, sg)) continue;

      // ignore zeros that come from cells outside graph convex hull
      if (r > 1e-4 && mu < 1e-12 && sg < 1e-12) continue;

      sum_mu += mu;
      sum_sg += sg;
      cnt++;
    }

    if (cnt > 0) {
      last_mu = sum_mu / cnt;
      last_sg = sum_sg / cnt;
    }

    m_mu_r_ns[ir]  = static_cast<float>(last_mu);
    m_sig_r_ns[ir] = static_cast<float>(last_sg);
  }
}


// 1D radial interpolation for r inside the table,
// and XTREL/GMCT-style extrapolation outside:
// - quadratic using last 3 points (your request)
// - if not enough points, linear using last 2 points
bool DCHXT2DLUT::radialMeanSigma(double r_cm, double& mu_ns, double& sig_ns) const
{
  if (m_rbin_cm.empty()) return false;

  const int n = (int)m_rbin_cm.size();
  if (n < 3) return false;

  if (r_cm <= m_rbin_cm.front()) {
    mu_ns  = m_mu_r_ns.front();
    sig_ns = m_sig_r_ns.front();
    return true;
  }

  // inside range: linear interpolation
  if (r_cm <= m_rbin_cm.back()) {
    int hi = std::lower_bound(m_rbin_cm.begin(), m_rbin_cm.end(), (float)r_cm) - m_rbin_cm.begin();
    if (hi <= 0) hi = 1;
    if (hi >= n) hi = n - 1;
    const int lo = hi - 1;

    const double r1 = m_rbin_cm[lo];
    const double r2 = m_rbin_cm[hi];
    const double t  = (r_cm - r1) / (r2 - r1);

    mu_ns  = (1.0 - t) * m_mu_r_ns[lo]  + t * m_mu_r_ns[hi];
    sig_ns = (1.0 - t) * m_sig_r_ns[lo] + t * m_sig_r_ns[hi];

    if (mu_ns < 0.0) mu_ns = 0.0;
    if (sig_ns < 0.0) sig_ns = 0.0;
    return true;
  }

  // outside range: extrapolation (like XTREL "after last point")
  const int i3 = n - 1;
  const int i2 = n - 2;
  const int i1 = n - 3;

  const double r1 = m_rbin_cm[i1];
  const double r2 = m_rbin_cm[i2];
  const double r3 = m_rbin_cm[i3];

  const double mu1 = m_mu_r_ns[i1];
  const double mu2 = m_mu_r_ns[i2];
  const double mu3 = m_mu_r_ns[i3];

  const double sg1 = m_sig_r_ns[i1];
  const double sg2 = m_sig_r_ns[i2];
  const double sg3 = m_sig_r_ns[i3];

  // quadratic extrapolation (3 last points)
  double mu_q = quad3(r_cm, r1, r2, r3, mu1, mu2, mu3);
  double sg_q = quad3(r_cm, r1, r2, r3, sg1, sg2, sg3);

  // linear fallback (2 last points) - same idea as XTREL
  const double slope_mu = (mu3 - mu2) / (r3 - r2);
  const double slope_sg = (sg3 - sg2) / (r3 - r2);
  const double mu_lin = mu3 + slope_mu * (r_cm - r3);
  const double sg_lin = sg3 + slope_sg * (r_cm - r3);

  mu_ns  = (std::isfinite(mu_q) ? mu_q : mu_lin);
  sig_ns = (std::isfinite(sg_q) ? sg_q : sg_lin);

  // enforce physical sanity:
  if (!std::isfinite(mu_ns) || mu_ns < 0.0) mu_ns = std::max(0.0, mu_lin);
  if (!std::isfinite(sig_ns) || sig_ns < 0.0) sig_ns = std::max(0.0, sg_lin);

  // do not let time decrease beyond last point
  if (mu_ns < mu3) mu_ns = std::max(mu3, mu_lin);

  return true;
}


bool DCHXT2DLUT::meanSigma(double x_cm, double y_cm, double& mu_ns, double& sig_ns) const
{
  if (!m_loaded) return false;

  // If inside the grid rectangle, use 2D bilinear value.
  // If outside, use 1D radial (with extrapolation beyond last point).
  if (bilinearAt(x_cm, y_cm, mu_ns, sig_ns)) {
    return true;
  }

  const double r = std::hypot(x_cm, y_cm);
  return radialMeanSigma(r, mu_ns, sig_ns);
}


double DCHXT2DLUT::sampleTimeNs(double x_cm, double y_cm, TRandom3& rng) const
{
  double mu = 0.0, sg = 0.0;
  if (!meanSigma(x_cm, y_cm, mu, sg)) return 0.0;

  double t = (sg > 0.0) ? rng.Gaus(mu, sg) : mu;
  if (t < 0.0) t = 0.0;
  return t;
}

