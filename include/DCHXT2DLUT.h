// @author Muhammad Saiel

#pragma once

#include <cmath>
#include <memory>
#include <string>
#include "TFile.h"

class TFile;
class TGraph2D;
class TRandom3;

class DCHXT2DLUT {
public:
  DCHXT2DLUT(const std::string& rootFile,
             const std::string& meanName = "xt_mean",
             const std::string& sigmaName = "xt_error");
  ~DCHXT2DLUT() = default;

  bool meanSigma(double x_cm, double y_cm, double& mu_ns, double& sig_ns) const;
  double sampleTimeNs(double x_cm, double y_cm, TRandom3& rng) const;

private:
  std::unique_ptr<TFile> m_file;
  TGraph2D* m_gMean{nullptr};
  TGraph2D* m_gSigma{nullptr};

  static bool isInvalid(double mu, double sig, double x, double y);
};

