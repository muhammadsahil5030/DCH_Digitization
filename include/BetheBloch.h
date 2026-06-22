#pragma once

#include <cmath>

namespace BB {
	inline double bethe_bloch(double p_GeV, double m_MeV, double I_MeV) {
  
		// Material / particle constants for He-based gas approximation
		constexpr double Z = 2.0;         // Atomic number of material
		constexpr double A = 4.0026;      // Relative atomic mass
		constexpr double z_particle = 1.0;
		constexpr double me_MeV = 0.511;  // Electron rest mass [MeV]
		constexpr double e = 1.0;
  		constexpr double K = 0.307075;
 
		// Density-effect parameters
	      	constexpr double c = 11.1393;
	      	constexpr double x0 = 2.2017;
	      	constexpr double x1 = 3.6122;
	      	constexpr double a = 0.13443;
	      	constexpr double k1 = 5.8347;
	      	constexpr double k2 = 2.0 * std::log(10.0);

	      	// Convert momentum from GeV to MeV
	      	const double p_MeV = p_GeV * 1.0e3;

	      	const double x = std::log10(p_MeV / m_MeV);
	      	double delta = 0.0;

	      	const double beta = p_MeV / std::sqrt(p_MeV * p_MeV + m_MeV * m_MeV);
	      	const double gamma = 1.0 / std::sqrt(1.0 - beta * beta);

	      	if (x >= x1) {
		    	delta = k2 * x - c;
	      	} else if (x >= x0) {
		    	delta = k2 * x - c + a * std::pow(x1 - x, k1);
	      	}

	      	// Maximum kinetic energy transferable to a free electron
	      	const double Tmax =
		  	2.0 * me_MeV * std::pow(gamma * beta, 2) /
		  	(1.0 + (2.0 * gamma * me_MeV / m_MeV) + std::pow(me_MeV / m_MeV, 2));

	      	// Bethe-Bloch formula
	      	const double dEdx =
		  	K * std::pow(z_particle * e, 2) * (Z / A) / std::pow(beta, 2) *
		  	(0.5 * std::log(2.0 * me_MeV * std::pow(beta * gamma, 2) * Tmax / (I_MeV * I_MeV)) -
		  	 std::pow(beta, 2) - 0.5 * delta);

	      	return dEdx;
	}
}  // namespace BB
