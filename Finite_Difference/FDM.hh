#ifndef FDM_HH
#define FDM_HH

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
//#include <Eigen/superlu>

#include "PDE.hh"

typedef std::vector<Eigen::VectorXd> matrix;

class FDMBase {
 protected:
  BlackScholesPDE* pde;

  //Space domain [Nminus*dx,Nplus*dx]
  int long Nminus = -500;
  unsigned long Nplus =500;
  unsigned long n;
  double dx;
  std::vector<double> x_values;

  //Time domain [0, 1/2*simga^2*T]
  unsigned long M; //number of time intervals
  double dt;

  //coefficients
  double alpha;

  //initial vector;
  Eigen::VectorXd u0;

  // Constructor
  FDMBase(unsigned long _n , unsigned long _M, BlackScholesPDE* _pde): n(_n), M(_M), pde(_pde){}

  virtual void calculate_step_sizes() = 0;
  virtual void set_initial_conditions() = 0;
  virtual void calculate_boundary_conditions_call(Eigen::VectorXd&) = 0;
  virtual void calculate_boundary_conditions_put() = 0;

 public:

  virtual matrix solve() = 0;
};


class FDMEulerExplicit : public FDMBase {
 protected:
  void calculate_step_sizes();
  void set_initial_conditions();
  void calculate_boundary_conditions_call(Eigen::VectorXd&);
  virtual void calculate_boundary_conditions_put();

 public:
  FDMEulerExplicit(unsigned long _n, unsigned long _M, BlackScholesPDE* _pde);

  matrix solve();
};

#endif
