#ifndef OPTION_HH
#define OPTION_HH

#include <iostream>
#include <algorithm>
#include <math.h>

class Option {
  protected:
    double E;
    double r;
    double T;
    double sigma;
    double k;
    double pay_off;


  public:
    Option() = default;

    virtual void set_payoff(double) = 0;

    //getter methods
    double
    get_E() const {return E;}

    double
    get_r() const {return r;}

    double
    get_T() const {return T;}

    double
    get_sigma() const {return sigma;}

    double
    get_payoff() const{return pay_off;}

    double
    get_k() const{return k;}


};


class EuropeanCall: public Option{
  public:
    EuropeanCall(double _E, double _r, double _T, double _sigma);

    void set_payoff(double x);

};


class EuropeanPut: public Option{
  public:
    EuropeanPut(double _E, double _r, double _T, double _sigma);

    void set_payoff(double x);

};

#endif
