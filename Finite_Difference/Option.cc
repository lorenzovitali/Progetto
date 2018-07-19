#include "Option.hh"

EuropeanCall::EuropeanCall(double _E, double _r, double _T, double _sigma){
  E = _E;
  r = _r;
  T = _T;
  sigma = _sigma;

  k = r/(0.5*sigma*sigma);
}

EuropeanPut::EuropeanPut(double _E, double _r, double _T, double _sigma){
  E = _E;
  r = _r;
  T = _T;
  sigma = _sigma;

  k = r/(0.5*sigma*sigma);
}

void EuropeanCall::set_payoff(double x){
  pay_off = std::max(exp(0.5*(k+1)*x) - exp(0.5*(k-1)*x), 0.0);
}

void EuropeanPut::set_payoff(double x){
  pay_off = std::max(exp(0.5*(k-1)*x) - exp(0.5*(k+1)*x), 0.0);
}
