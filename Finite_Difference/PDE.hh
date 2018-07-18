#ifndef PDE_HH
#define PDE_HH

#include "Option.hh"

class BlackScholesPDE {
  private:
    Option* option;

  public:
    BlackScholesPDE(Option* _option);

    //coefficients
    double diff_coeff(void) const;

    //boundaries for CALL
    double
    call_boundary_left(void) const;

    double
    call_boundary_right(double S_max) const;

    double
    put_boundary_left(double t) const;

    double
    put_boundary_right(void) const;

    double
    init_cond(double x);




};

#endif
