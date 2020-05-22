// Copyright 2020 Carmelo Evoli (GSSI) - MIT License
#ifndef INCLUDE_BREMSSTRAHLUNG_FUNCTIONS_H_
#define INCLUDE_BREMSSTRAHLUNG_FUNCTIONS_H_

namespace GALPROP {

#define pow2(A) ((A) * (A))
#define pow3(A) ((A) * (A) * (A))
#define pow4(A) ((A) * (A) * (A) * (A))

double Elwert_factor(double beta_i, double beta_f, int Z);
double xi_func(double T_electron, double k, int Z, int N);
double Phi_u(double gamma_i, double gamma_f, double k);
double Phi_1(double gamma_i, double gamma_f, double k, double delta, int Z, int N);
double Phi_2(double gamma_i, double gamma_f, double k, double delta, int Z, int N);

}  // namespace GALPROP

#endif  // INCLUDE_BREMSSTRAHLUNG_FUNCTIONS_H_
