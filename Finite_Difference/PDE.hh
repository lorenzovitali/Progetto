#ifndef PDE_HH
#define PDE_HH

#include <memory>
#include "Option.hh"

class BlackScholesPDE {
  private:
    //std::shared_ptr<Option> option;
    Option* option;

  public:
    BlackScholesPDE(Option* _option): option(_option){}


    //coefficients
    double diff_coeff(void) const;

    Option* get_option()const{return option;}

    //boundaries for CALL
    double
    call_boundary_left(void) const;

    double
    call_boundary_right(double x_max, double tau) const;

    double
    put_boundary_left(double tau) const;

    double
    put_boundary_right(void) const;

    double
    init_cond(double x);




};

#endif
