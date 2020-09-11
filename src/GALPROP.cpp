// Copyright 2020 Carmelo Evoli (GSSI) - MIT License
#include "include/GALPROP.h"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>

#include <cassert>
#include <cmath>

#include "include/cgs.h"

#define LIMIT 10000
#define EPSINT 1e-5
#define KEYINT 3

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
    double value = 16. * pow2(static_cast<double>(Z) * cgs::electron_radius) * cgs::alpha_f;
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

    double value = pow2(static_cast<double>(Z) * cgs::electron_radius) * cgs::alpha_f *
                   (p_f / p_i) / k * factor;
    return value;
}

double BremsstrahlungSpectrum::dsdk_high_energy() {
    double phi_1 =
        (N == 0) ? pow2(Z) * Phi_u(gamma_i, gamma_f, k) : Phi_1(gamma_i, gamma_f, k, delta, Z, N);
    double phi_2 =
        (N == 0) ? pow2(Z) * Phi_u(gamma_i, gamma_f, k) : Phi_2(gamma_i, gamma_f, k, delta, Z, N);
    double factor = (1. + pow2(gamma_f / gamma_i)) * phi_1 - 2. / 3. * gamma_f / gamma_i * phi_2;
    double value = pow2(cgs::electron_radius) * cgs::alpha_f / k * factor;
    return value;
}

double Elwert_factor(double beta_i, double beta_f, int Z) {
    double value =
        beta_i * (1. - std::exp(-2. * M_PI * static_cast<double>(Z) * cgs::alpha_f / beta_i));
    value /= beta_f * (1. - std::exp(-2. * M_PI * static_cast<double>(Z) * cgs::alpha_f / beta_f));
    return value;
}

double xi_func(double T_electron, double k, int Z, int N) {
    double b = 0.07 * cgs::MeV;
    double c = 0.33 * cgs::MeV;
    double xi = 1. + static_cast<double>(N) / pow2(Z) * (1. - std::exp((b - T_electron) / 9. / b));
    xi *= 1. - 0.3 * std::exp(-k / c);
    return xi;
}

struct gslParams_I_Phi_N {
    double delta;
    int Z;
    int N;
};

double R_1(double q, double Z) {
    double F_1 = 1. / pow2(1. + pow2(q) / pow2(2. * cgs::alpha_f * static_cast<double>(Z)));
    return 1 - F_1;
}

double R_2(double q, double Z) {
    double F_2 =
        1. / pow2(1. + pow2(q) / pow2(2. * cgs::alpha_f * (static_cast<double>(Z) - 5. / 16.)));
    return 2. * (1. - F_2) - (1. - pow2(F_2)) / static_cast<double>(Z);
}

double gslFunc_I_Phi_1(double q, void *params) {
    gslParams_I_Phi_N prms = *(reinterpret_cast<gslParams_I_Phi_N *>(params));
    double R_N = (prms.N == 1) ? R_1(q, prms.Z) : R_2(q, prms.Z);
    return R_N / pow3(q) * pow2(q - prms.delta);
}

double I_Phi_1(double delta, int Z, int N) {
    double result, error;
    gslParams_I_Phi_N params = {delta, Z, N};

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(LIMIT);
    gsl_function F;
    F.function = &gslFunc_I_Phi_1;
    F.params = &params;
    gsl_integration_qags(&F, delta, 1, 0, EPSINT, LIMIT, w, &result, &error);
    gsl_integration_workspace_free(w);

    return result;
}

double gslFunc_I_Phi_2(double q, void *params) {
    gslParams_I_Phi_N prms = *(reinterpret_cast<gslParams_I_Phi_N *>(params));
    double R_N = (prms.N == 1) ? R_1(q, prms.Z) : R_2(q, prms.Z);
    double factor = pow3(q) - 6. * pow2(prms.delta) * q * std::log(q / prms.delta) +
                    3. * pow2(prms.delta) * q - 4. * pow3(prms.delta);
    double value = R_N / pow4(q) * factor;
    return value;
}

double I_Phi_2(double delta, int Z, int N) {
    double result, error;
    gslParams_I_Phi_N params = {delta, Z, N};

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(LIMIT);
    gsl_function F;
    F.function = &gslFunc_I_Phi_2;
    F.params = &params;
    gsl_integration_qags(&F, delta, 1, 0, EPSINT, LIMIT, w, &result, &error);
    gsl_integration_workspace_free(w);

    return result;
}

double Phi_u(double gamma_i, double gamma_f, double k) {
    return 4. * (std::log(2. * gamma_i * gamma_f / k) - 0.5);
}

double Phi_1(double gamma_i, double gamma_f, double k, double delta, int Z, int N) {
    double I = I_Phi_1(delta, Z, N);
    double value = static_cast<double>(pow2(Z - N)) * Phi_u(gamma_i, gamma_f, k) +
                   8. * Z * (1. - (N - 1.) / static_cast<double>(Z) + I);
    return value;
}

double Phi_2(double gamma_i, double gamma_f, double k, double delta, int Z, int N) {
    double I = I_Phi_2(delta, Z, N);
    double value = static_cast<double>(pow2(Z - N)) * Phi_u(gamma_i, gamma_f, k) +
                   8. * Z * (5. / 6. * (1. - (N - 1.) / static_cast<double>(Z)) + I);
    return value;
}

}  // namespace GALPROP
