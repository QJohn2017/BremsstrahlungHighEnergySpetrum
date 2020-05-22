// Copyright 2020 Carmelo Evoli (GSSI) - MIT License
#include "include/bremsstrahlung_sigma.h"

#include <cassert>
#include <cmath>

#include "include/bremsstrahlung_functions.h"
#include "include/cgs.h"

namespace GALPROP {

BremsstrahlungSpectrum::BremsstrahlungSpectrum(int atomic_number, int atomic_electrons)
    : Z(atomic_number), N(atomic_electrons) {
    assert(atomic_number == 1 || atomic_number == 2);
    assert(atomic_electrons == 1 || atomic_electrons == 2);
}

void BremsstrahlungSpectrum::set_E_gamma(double E_gamma) {
    this->E_gamma = E_gamma;
    k = E_gamma / cgs::electron_mass_c2;
}

void BremsstrahlungSpectrum::set_T_electron(double T_electron) {
    this->T_electron = T_electron;
    double T_electron_final = T_electron - E_gamma;
    p_i = std::sqrt(T_electron * (T_electron + 2. * cgs::electron_mass_c2)) / cgs::electron_mass_c2;
    p_f = std::sqrt(T_electron_final * (T_electron_final + 2. * cgs::electron_mass_c2)) /
          cgs::electron_mass_c2;
    gamma_i = (T_electron + cgs::electron_mass_c2) / cgs::electron_mass_c2;
    gamma_f = (T_electron_final + cgs::electron_mass_c2) / cgs::electron_mass_c2;
    beta_i = std::sqrt(1. - 1. / pow2(gamma_i));
    beta_f = std::sqrt(1. - 1. / pow2(gamma_f));
}

void BremsstrahlungSpectrum::init_kinematic(double T_electron, double E_gamma) {
    set_E_gamma(E_gamma);
    set_T_electron(T_electron);
    delta = k / 2. / gamma_i / gamma_f;
}

double BremsstrahlungSpectrum::get(double T_electron, double E_gamma) {
    init_kinematic(T_electron, E_gamma);

    if (T_electron > E_gamma) {
        if (T_electron < 0.07 * cgs::MeV) {
            return Elwert_factor(beta_i, beta_f, Z) * dsdk_low_energy();
        } else if (T_electron < 2. * cgs::MeV) {
            return xi_func(T_electron, k, Z, N) * Elwert_factor(beta_i, beta_f, Z) *
                   dsdk_intermediate_energy();
        } else {
            return dsdk_high_energy();
        }
    } else {
        return 0;
    }
}

double BremsstrahlungSpectrum::dsdk_low_energy() {
    double value = 16. * pow2(static_cast<double>(Z) * cgs::r_e) * cgs::alpha_f;
    value /= 3. * k * pow2(p_i);
    value *= std::log((p_i + p_f) / (p_i - p_f));
    return value;
}

double BremsstrahlungSpectrum::dsdk_intermediate_energy() {
    double L = 2. * std::log((gamma_i * gamma_f + p_i * p_f - 1.) / k);
    double epsilon_i = std::log((gamma_i + p_i) / (gamma_i - p_i));
    double epsilon_f = std::log((gamma_f + p_f) / (gamma_f - p_f));

    double ininner_factor = epsilon_i * (gamma_i * gamma_f + pow2(p_i)) / pow3(p_i);
    ininner_factor -= epsilon_f * (gamma_i * gamma_f + pow2(p_f)) / pow3(p_f);
    ininner_factor += 2. * k * gamma_i * gamma_f / pow2(p_i * p_f);

    double inner_factor = 8. / 3. * gamma_i * gamma_f / p_i / p_f;
    inner_factor += pow2(k) * (pow2(gamma_i * gamma_f) + pow2(p_i * p_f)) / pow3(p_i * p_f);
    inner_factor += k / 2. / p_i / p_f * ininner_factor;

    double factor = 4. / 3.;
    factor -= 2. * gamma_i * gamma_f * (pow2(p_f) + pow2(p_i)) / pow2(p_f * p_i);
    factor += epsilon_i * gamma_f / pow3(p_i) + epsilon_f * gamma_i / pow3(p_f);
    factor -= epsilon_i * epsilon_f / p_i / p_f;
    factor += L * inner_factor;

    double value =
        pow2(static_cast<double>(Z) * cgs::r_e) * cgs::alpha_f * (p_f / p_i) / k * factor;
    return value;
}

double BremsstrahlungSpectrum::dsdk_high_energy() {
    double phi_1 =
        (N == 0) ? pow2(Z) * Phi_u(gamma_i, gamma_f, k) : Phi_1(gamma_i, gamma_f, k, delta, Z, N);
    double phi_2 =
        (N == 0) ? pow2(Z) * Phi_u(gamma_i, gamma_f, k) : Phi_2(gamma_i, gamma_f, k, delta, Z, N);
    double factor = (1. + pow2(gamma_f / gamma_i)) * phi_1 - 2. / 3. * gamma_f / gamma_i * phi_2;
    double value = pow2(cgs::r_e) * cgs::alpha_f / k * factor;
    return value;
}

}  // namespace GALPROP
