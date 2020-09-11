// Copyright 2020 Carmelo Evoli (GSSI) - MIT License
#ifndef INCLUDE_LOSSRATE_H_
#define INCLUDE_LOSSRATE_H_

#include <include/cgs.h>

#include <cassert>
#include <cmath>

class LossRate {
   public:
    LossRate() {}
    virtual ~LossRate() = default;
    double dEdt_neutral_gas(double T, double n_HI) const;
    double dEdt_ionized_gas(double T, double Z, double n_HII) const;

    double get(double T) const {
        return dEdt_neutral_gas(T, 1. / cgs::cm3) + dEdt_ionized_gas(T, 1, 0.01);
    }

   private:
    const double f_He = 1.4;
    const double T_H = 62.8 * cgs::gram / cgs::cm2;
    const double T_He = 93.1 * cgs::gram / cgs::cm2;
};

#endif  // INCLUDE_LOSSRATE_H_
