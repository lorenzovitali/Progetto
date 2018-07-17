#ifndef OPTION_HH
#define OPTION_HH

#include "payoff.h"

class Option {
 public:
  PayOff* pay_off;

  double K;
  double r;
  double T;
  double sigma;

  VanillaOption();
  VanillaOption(double _K, double _r, double _T,
                double _sigma, PayOff* _pay_off);
};

#endif
