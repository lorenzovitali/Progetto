#ifndef FDM_HH
#define FDM_HH

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include "PDE.hh"

typedef std::vector<Eigen::VectorXd> matrix;

class FDMBase {
 protected:
  BlackScholesPDE* pde;

  //Space domain [Nminus*dx,Nplus*dx]
  int long Nminus = -500;
  unsigned long Nplus =500;
  unsigned long n; //number of space intervals
  double dx;
  std::vector<double> x_values;

  //Time domain [0, 1/2*simga^2*T]
  unsigned long M; //number of time intervals
  double dt;

  //coefficients
  double alpha;

  //initial vector;
  Eigen::VectorXd u0;

  //matrices
  Eigen::SparseMatrix<double, Eigen::RowMajor> A, I;

  // Constructor
  FDMBase(unsigned long _n , unsigned long _M, BlackScholesPDE* _pde);

  virtual void calculate_step_sizes();
  virtual void set_initial_conditions();
  virtual void calculate_boundary_conditions_call(Eigen::VectorXd&);
  //virtual void calculate_boundary_conditions_put();

 public:

  void BuildMatrix();
  virtual matrix solve() = 0;
};


class FDMEulerExplicit : public FDMBase {
 public:
  FDMEulerExplicit(unsigned long _n, unsigned long _M, BlackScholesPDE* _pde):FDMBase(_n, _M, _pde){}

  matrix solve();
};


class FDMEulerImplicit : public FDMBase{
  public:
    FDMEulerImplicit(unsigned long _n, unsigned long _M, BlackScholesPDE* _pde):FDMBase(_n, _M, _pde){}

    matrix solve();
};

class FDMCranckNicholson : public FDMBase{
public:
  FDMCranckNicholson(unsigned long _n, unsigned long _M, BlackScholesPDE* _pde):FDMBase(_n, _M, _pde){}

  matrix solve();
};




#endif
