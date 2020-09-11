// Copyright 2020 Carmelo Evoli (GSSI) - MIT License
#ifndef INCLUDE_BREMSSTRAHLUNG_SIGMA_H_
#define INCLUDE_BREMSSTRAHLUNG_SIGMA_H_

namespace GALPROP {

double Elwert_factor(double beta_i, double beta_f, int Z);
double xi_func(double T_electron, double k, int Z, int N);
double Phi_u(double gamma_i, double gamma_f, double k);
double Phi_1(double gamma_i, double gamma_f, double k, double delta, int Z, int N);
double Phi_2(double gamma_i, double gamma_f, double k, double delta, int Z, int N);

class BremsstrahlungSpectrum {
   public:
    BremsstrahlungSpectrum(int atomic_number, int atomic_electrons);
    virtual ~BremsstrahlungSpectrum() {}
    void set_E_gamma(double E_gamma);
    void set_T_electron(double T_electron);
    void init_kinematic(double T_electron, double E_gamma);
    double get(double T_electron, double E_gamma);

   protected:
    int Z;
    int N;

    double E_gamma = -1;
    double k = -1;
    double T_electron = -1;
    double p_i = -1;
    double p_f = -1;
    double beta_i = -1;
    double beta_f = -1;
    double gamma_i = -1;
    double gamma_f = -1;
    double delta = -1;

   private:
    double dsdk_low_energy();
    double dsdk_intermediate_energy();
    double dsdk_high_energy();
};

}  // namespace GALPROP

#endif  // INCLUDE_BREMSSTRAHLUNG_SIGMA_H_
