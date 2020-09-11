// Copyright 2020 Carmelo Evoli (GSSI) - MIT License
#include "include/lossrate.h"

/*
 * ref: Strong & Moskalenko Apj, 1998, (C12)
 */
double LossRate::dEdt_neutral_gas(double T, double n_HI) const {
    const double E = T + cgs::electron_mass_c2;
    const double gamma_e = E / cgs::electron_mass_c2;
    assert(gamma_e > 800.);
    double value = n_HI * cgs::proton_mass * (1. / T_H + 4. * f_He / T_He);
    value *= cgs::c_light * E;
    return value;
}

/*
 * ref: Strong & Moskalenko Apj, 1998, (C10)
 */
double LossRate::dEdt_ionized_gas(double T, double Z, double n_HII) const {
    const double E = T + cgs::electron_mass_c2;
    const double gamma_e = E / cgs::electron_mass_c2;
    constexpr double factor = cgs::alpha_f * pow2(cgs::electron_radius) * cgs::c_light;
    double value = 4. * factor * (Z * (Z + 1)) * n_HII * E;
    value *= std::log(2. * gamma_e) - 1. / 3.;
    return value;
}
