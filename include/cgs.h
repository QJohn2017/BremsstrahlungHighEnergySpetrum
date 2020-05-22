// Copyright 2020 Carmelo Evoli (GSSI) - MIT License
#ifndef INCLUDE_CGS_H_
#define INCLUDE_CGS_H_

namespace cgs {

// base
const double cm = 1;
const double gr = 1;
const double erg = 1;
const double sec = 1;

// derived
const double cm2 = cm * cm;
const double eV = 1.6022e-12 * erg;
const double MeV = 1e6 * eV;
const double GeV = 1e3 * MeV;
const double TeV = 1e3 * GeV;
const double r_e = 2.81794e-13 * cm;
const double alpha_f = 1. / 137.;
const double electron_mass_c2 = 0.511 * MeV;
const double mbarn = 1e-27 * cm * cm;
const double h = 6.626e-27 * cm2 * gr / sec;
const double c = 2.99792458e10 * cm / sec;
const double hc = h * c;

}  // namespace cgs

#endif  // INCLUDE_CGS_H_
