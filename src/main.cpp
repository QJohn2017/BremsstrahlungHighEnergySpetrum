// Copyright 2020 Carmelo Evoli (GSSI) - MIT License
#include <iomanip>
#include <iostream>

#include "include/GALPROP.h"
#include "include/Tsai74.h"
#include "include/cgs.h"
#include "include/lossrate.h"

void print_spectrum() {
    GALPROP::BremsstrahlungSpectrum galpropSpectrum(1, 1);
    TSAI74::Tsai74 hermesSpectrum(1);
    std::cout << std::scientific;
    for (double E_gamma = 0.01 * cgs::MeV; E_gamma < 1e4 * cgs::MeV; E_gamma *= 1.1) {
        std::cout << E_gamma / cgs::MeV << "\t";
        double k = E_gamma / cgs::electron_mass_c2;
        std::cout << k * galpropSpectrum.get(1e3 * cgs::MeV, E_gamma) / cgs::mbarn << "\t";
        std::cout << hermesSpectrum.get(1e3 * cgs::MeV, E_gamma) / cgs::mbarn << "\t";
        std::cout << "\n";
    }
}

double CMB_timescale(double T, double energy_density) {
    const double gamma_electron = (T + cgs::electron_mass_c2) / cgs::electron_mass_c2;
    const double factor = 3. * cgs::electron_mass_c2 / 4. / cgs::sigma_th / cgs::c_light;
    return factor / energy_density / gamma_electron;
}

void print_lossrate() {
    LossRate dEdt;
    for (double T = cgs::GeV; T < 1e1 * cgs::GeV; T *= 1.01) {
        std::cout << T / cgs::GeV << " ";
        std::cout << dEdt.get(T) / (cgs::GeV / cgs::sec) << " ";
        std::cout << 50. * T / dEdt.get(T) / (cgs::year) << " ";
        std::cout << CMB_timescale(T, 0.25 * cgs::eV / cgs::cm3) / (cgs::year) << " ";
        std::cout << "\n";
    }
}

int main() { print_lossrate(); }
