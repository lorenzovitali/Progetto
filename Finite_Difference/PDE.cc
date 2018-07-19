#include "PDE.hh"

// Diffusion coefficient
double BlackScholesPDE::diff_coeff(void) const {
  return 1;
}

double BlackScholesPDE::call_boundary_left()const{
  return 0.0;
}

double BlackScholesPDE::call_boundary_right(double S_max) const {
  return S_max;
}

double BlackScholesPDE::put_boundary_left(double t) const{
  return option->get_E()*exp(-option->get_r() * (option->get_T() - t));
}

double BlackScholesPDE::put_boundary_right(void)const{
  return 0.0;
}

// Initial condition
double BlackScholesPDE::init_cond(double x){
  option->set_payoff(x);
  return option->get_payoff();
}
