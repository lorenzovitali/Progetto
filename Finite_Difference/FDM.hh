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
  std::shared_ptr<BlackScholesPDE> pde;

  //Space domain [Nminus*dx,Nplus*dx]
  int long Nminus;
  int long Nplus;
  int long N; //number of space intervals
  double dx;
  std::vector<double> x_values;

  //Time domain [0, 1/2*simga^2*T]
  unsigned long M; //number of time intervals
  double dt;
  std::vector<double> tau_values;

  //coefficients
  double alpha;

  //initial vector;
  Eigen::VectorXd u0;

  //matrices
  Eigen::SparseMatrix<double, Eigen::RowMajor> A, I;

  // Constructor
  FDMBase(int _N , double dx, unsigned long _M, std::shared_ptr<BlackScholesPDE> _pde);

  void calculate_step_sizes();
  void set_initial_conditions();
  void calculate_boundary_conditions(Eigen::VectorXd&, double);

 public:

  double get_dt()const{ return dt; }
  double get_dx()const{ return dx; }
  double get_alpha()const{ return alpha; }
  std::vector<double> get_x() const{ return x_values; }
  std::vector<double> get_tau() const{ return tau_values; }


  void BuildMatrix();
  void change_var(matrix&);
  virtual matrix solve() = 0;
};


class FDMEulerExplicit : public FDMBase {
 public:
  FDMEulerExplicit(int _N, double dx, unsigned long _M, std::shared_ptr<BlackScholesPDE> _pde):
    FDMBase(_N, dx, _M, _pde){}

  matrix solve();
};


class FDMEulerImplicit : public FDMBase{
  public:
    FDMEulerImplicit(int _N, double dx, unsigned long _M, std::shared_ptr<BlackScholesPDE> _pde):
      FDMBase(_N, dx, _M, _pde){}

    matrix solve();
};

class FDMCranckNicholson : public FDMBase{
public:
  FDMCranckNicholson(int _N, double dx, unsigned long _M, std::shared_ptr<BlackScholesPDE> _pde):
    FDMBase(_N, dx, _M, _pde){}

  matrix solve();
};




#endif
