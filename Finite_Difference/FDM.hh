#ifndef FDM_HH
#define FDM_HH

#include "PDE.hh"
#include <vector>

class FDMBase {
 protected:
  BlackScholesPDE* pde;

  //Space domain [-N*dx,N*dx]
  unsigned long N = 1000;
  unsigned long n;
  double dx;
  std::vector<double> x_values;

  //Time domain [0, 1/2*simga^2*T]
  unsigned long M; //number of time intervals
  double dt;
  double prev_t, cur_t;

  //coefficients
  double alpha;

  // Constructor
  FDMBase(unsigned long _n , unsigned long _M, BlackScholesPDE* _pde);

  virtual void calculate_step_sizes() = 0;
  virtual void set_initial_conditions() = 0;
  virtual void calculate_boundary_conditions() = 0;
  virtual void calculate_inner_domain() = 0;

 public:
  // Carry out the actual time-stepping
  virtual void step_march() = 0;
};


class FDMEulerExplicit : public FDMBase {
 protected:
  void calculate_step_sizes();
  void set_initial_conditions();
  void calculate_boundary_conditions();
  void calculate_inner_domain();

 public:
  FDMEulerExplicit(unsigned long _n, unsigned long _M, BlackScholesPDE* _pde);

  void step_march();
};

#endif
