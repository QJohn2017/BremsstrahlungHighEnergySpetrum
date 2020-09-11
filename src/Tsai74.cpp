// Copyright 2020 Carmelo Evoli (GSSI) - MIT License
#include "include/Tsai74.h"

#include <gsl/gsl_integration.h>

#include <cmath>
#include <iostream>
#include <utility>

#include "include/cgs.h"

#define pow2(A) ((A) * (A))

namespace TSAI74 {

std::pair<double, double> RadiationLogarithms(int Z) {
    // elastic and inelatic radiation logarithms for light elements (where the
    // Thomas-Fermi model doesn't work): computed by using Dirac-Fock model of atom.
    const double gFelLowZet[] = {0.0, 5.3104, 4.7935, 4.7402, 4.7112, 4.6694, 4.6134, 4.5520};
    const double gFinelLowZet[] = {0.0, 5.9173, 5.6125, 5.5377, 5.4728, 5.4174, 5.3688, 5.3236};

    double Fel, Finel;
    if (Z < 5) {
        Fel = gFelLowZet[Z];
        Finel = gFinelLowZet[Z];
    } else {
        const double LogZ = std::log(static_cast<double>(Z));
        Fel = std::log(184.15) - LogZ / 3.;
        Finel = std::log(1194) - 2. * LogZ / 3.;
    }
    return {Fel, Finel};
}

double ComputeCoulombCorrection(int Z, size_t n_max = 100) {
    double value = 0;
    for (size_t n = 1; n < n_max; ++n)
        value += 1. / n / (n * n + cgs::alpha_f * cgs::alpha_f * static_cast<double>(Z * Z));
    return value;
}

Tsai74::Tsai74(int Z) : m_Z(Z) {}

double Tsai74::get(double primaryKineticEnergy, double gammaEnergy) {
    double xsecs = 0;

    if (gammaEnergy / primaryKineticEnergy < 1.) {
        xsecs = 4. * cgs::alpha_f * pow2(cgs::electron_radius) / 3.;

        const auto PrimaryTotalEnergy = primaryKineticEnergy + cgs::electron_mass_c2;
        const auto y = gammaEnergy / PrimaryTotalEnergy;
        const auto f = ComputeCoulombCorrection(m_Z);
        const auto radLogs = RadiationLogarithms(m_Z);
        const auto F_el = radLogs.first;
        const auto F_inel = radLogs.second;

        const auto X1 = (4. - 4. * y + y * y / 3.) * (m_Z * m_Z * (F_el - f) + m_Z * F_inel);
        const auto X2 = (1. - y) * (m_Z * m_Z + m_Z) / 3.;

        xsecs *= X1 + X2;
    }
    return xsecs;
}

}  // namespace TSAI74
