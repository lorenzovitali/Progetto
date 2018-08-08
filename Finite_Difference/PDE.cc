#include "PDE.hh"

double BlackScholesPDE::call_boundary_left()const{
  return 0.0;
}

double BlackScholesPDE::call_boundary_right(double x_max, double tau) const {
  return exp(0.5 * (option->get_k() + 1 ) * x_max + 0.25 * (option->get_k() + 1 )*(option->get_k() + 1 ) * tau);
}

double BlackScholesPDE::put_boundary_left(double x_min, double tau) const{
  return exp(0.5 * (option->get_k() - 1 ) * x_min + 0.25 * (option->get_k() - 1 )*(option->get_k() - 1 ) * tau);
}

double BlackScholesPDE::put_boundary_right(void)const{
  return 0.0;
}

// Initial condition
double BlackScholesPDE::init_cond(double x){
  option->set_payoff(x);
  return option->get_payoff();
}
