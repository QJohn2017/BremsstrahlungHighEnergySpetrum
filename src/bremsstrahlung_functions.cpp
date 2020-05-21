// Copyright 2020 Carmelo Evoli (GSSI) - MIT License
#include "include/bremsstrahlung_functions.h"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>

#include <cmath>

#include "include/cgs.h"

#define LIMIT 10000
#define EPSINT 1e-5
#define KEYINT 3

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
