// Copyright 2020 Carmelo Evoli (GSSI) - MIT License
#include <iomanip>
#include <iostream>

#include "include/bremsstrahlung_sigma.h"
#include "include/cgs.h"

int main() {
    BremsstrahlungSpectrum spectrum(1, 1);

    std::cout << std::scientific;
    for (double E_gamma = 0.01 * cgs::MeV; E_gamma < 1e4 * cgs::MeV; E_gamma *= 1.1) {
        std::cout << E_gamma / cgs::MeV << "\t";
        double k = E_gamma / cgs::electron_mass_c2;
        std::cout << k * spectrum.get(1e1 * cgs::MeV, E_gamma) / cgs::mbarn << "\t";
        std::cout << k * spectrum.get(1e2 * cgs::MeV, E_gamma) / cgs::mbarn << "\t";
        std::cout << k * spectrum.get(1e3 * cgs::GeV, E_gamma) / cgs::mbarn << "\t";
        std::cout << k * spectrum.get(1e4 * cgs::GeV, E_gamma) / cgs::mbarn << "\t";
        std::cout << "\n";
    }

    return 0;
}
