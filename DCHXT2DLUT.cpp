#include "DCHXT2DLUT.h"

#include "TFile.h"
#include "TGraph2D.h"

#include <iostream>
#include <limits>

#include "TRandom3.h"


namespace {

// Same “walk inward” idea you used, but ONLY during initialize when building the grid.
static bool queryGraphMeanSigma(TGraph2D& gMean,
                                TGraph2D& gSig,
                                double x_cm_in, double y_cm_in,
                                double& mu_ns, double& sig_ns)
{
  double x = x_cm_in;
  double y = y_cm_in;

  auto interp = [&](double xx, double yy, double& mu, double& sg) {
    mu = gMean.Interpolate(xx, yy);
    sg = gSig .Interpolate(xx, yy);
  };

  interp(x, y, mu_ns, sig_ns);

  auto invalid = [&](double mu, double sg, double xx, double yy) {
    if (!std::isfinite(mu) || !std::isfinite(sg)) return true;
    if (mu < 0.0 || sg < 0.0) return true;

    // If both are ~0 away from origin, treat as "out-of-domain"
    const double r = std::hypot(xx, yy);
    if (r > 1e-4 && std::abs(mu) < 1e-9 && std::abs(sg) < 1e-9) return true;

    return false;
  };

  double step = 1e-3; // cm
  for (int it = 0; it < 30 && invalid(mu_ns, sig_ns, x, y); ++it) {
    const double r = std::hypot(x, y);
    if (r < 1e-12) break;
    x -= step * x / r;
    y -= step * y / r;
    interp(x, y, mu_ns, sig_ns);
    step *= 2.0;
  }

  return (std::isfinite(mu_ns) && std::isfinite(sig_ns) && mu_ns >= 0.0 && sig_ns >= 0.0);
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

  if (nx < 10 || ny < 10) return false;

  TFile f(rootFile.c_str(), "READ");
  if (f.IsZombie()) return false;

  TGraph2D* gM = nullptr;
  TGraph2D* gS = nullptr;
  f.GetObject(meanName.c_str(),  gM);
  f.GetObject(sigmaName.c_str(), gS);
  if (!gM || !gS) return false;

  // Compute bounds from points
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

  // A tiny padding helps avoid edge issues
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

  // Build the grid NOW (single-thread, initialize)
  for (int iy = 0; iy < m_ny; ++iy) {
    const double y = m_ymin_cm + iy * m_dy;
    for (int ix = 0; ix < m_nx; ++ix) {
      const double x = m_xmin_cm + ix * m_dx;

      double mu = 0.0, sg = 0.0;
      const bool ok = queryGraphMeanSigma(*gM, *gS, x, y, mu, sg);

      if (!ok) {
        // Last-resort clamp slightly inward by rmax (still consistent with your older logic)
        const double r = std::hypot(x, y);
        if (r > 1e-12 && m_rmax_cm > 0.0) {
          const double scale = 0.999 * m_rmax_cm / r;
          queryGraphMeanSigma(*gM, *gS, x * scale, y * scale, mu, sg);
        } else {
          mu = 0.0; sg = 0.0;
        }
      }

      if (!std::isfinite(mu) || mu < 0.0) mu = 0.0;
      if (!std::isfinite(sg) || sg < 0.0) sg = 0.0;

      m_mu[idx(ix, iy)]  = static_cast<float>(mu);
      m_sig[idx(ix, iy)] = static_cast<float>(sg);
    }
  }

  m_loaded = true;
  return true;
}

bool DCHXT2DLUT::meanSigma(double x_cm, double y_cm, double& mu_ns, double& sig_ns) const
{
  if (!m_loaded) return false;

  // Clamp to grid bounds
  const double x = clamp(x_cm, m_xmin_cm, m_xmax_cm);
  const double y = clamp(y_cm, m_ymin_cm, m_ymax_cm);

  const double fx = (x - m_xmin_cm) / m_dx;
  const double fy = (y - m_ymin_cm) / m_dy;

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

  // Bilinear interpolation
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

double DCHXT2DLUT::sampleTimeNs(double x_cm, double y_cm, TRandom3& rng) const
{
  double mu = 0.0, sg = 0.0;
  if (!meanSigma(x_cm, y_cm, mu, sg)) return 0.0;

  double t = (sg > 0.0) ? rng.Gaus(mu, sg) : mu;
  if (t < 0.0) t = 0.0;
  return t;
}

