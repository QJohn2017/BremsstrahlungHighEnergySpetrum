// Copyright 2020 Carmelo Evoli (GSSI) - MIT License
#ifndef INCLUDE_TSAI74_H_
#define INCLUDE_TSAI74_H_

namespace TSAI74 {

class Tsai74 {
   public:
    explicit Tsai74(int Z);

    ~Tsai74() = default;

    double get(double T_electron, double E_gamma);

   protected:
    int m_Z = 1;
};

}  // namespace TSAI74

#endif  // INCLUDE_TSAI74_H_
