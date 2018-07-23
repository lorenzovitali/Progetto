#include <iostream>
#include "FDM.hh"

FDMBase::FDMBase(unsigned long _n, unsigned long _M, BlackScholesPDE* _pde):n(_n), M(_M), pde(_pde){

                                      calculate_step_sizes();
                                      set_initial_conditions();

                                    }

void FDMBase::BuildMatrix(){
  //bulding matrix A and I;
  A.resize(n,n);
  I.resize(n,n);

  using T = Eigen::Triplet<double>;
  std::vector<T> coeff;

  coeff.emplace_back(0, 0, 2);
  coeff.emplace_back(0, 1, -1);
  coeff.emplace_back(n-1 ,n-1, -1);
  coeff.emplace_back(n-1 ,n-2, 2);

  for (unsigned i = 1; i < n-1; i++){
    coeff.emplace_back(i,i,2);
    coeff.emplace_back(i,i-1,-1);
    coeff.emplace_back(i,i+1,-1);
  }

  A.setFromTriplets(coeff.begin(), coeff.end());
  coeff.clear();

  for (unsigned i = 0; i < n; i++){
    coeff.emplace_back(i,i,1);
  }

  I.setFromTriplets(coeff.begin(), coeff.end());
}


void FDMBase::calculate_step_sizes() {

  dx = (Nplus - Nminus) / n;

  dt = (0.5* pde->get_option()->get_sigma() * pde->get_option()->get_sigma() * pde->get_option()->get_T())/M;

}


void FDMBase::set_initial_conditions() {
  // Spatial settings
  double cur_spot = -Nminus*dx;
  u0.resize(M);

  for (unsigned long i=0; i<n; i++) {

    cur_spot = static_cast<double>(i)*dx;
    u0[i] = pde->init_cond(cur_spot);
    x_values.push_back(cur_spot);

  }
}


void FDMBase::calculate_boundary_conditions_call(Eigen::VectorXd & u_new) {
  u_new[0] = pde->call_boundary_left();
  u_new[n-1] = pde->call_boundary_right(Nplus*dx);
}


matrix FDMEulerExplicit::solve(){
  matrix result;
  BuildMatrix(); //building A and I;
  Eigen::VectorXd u = u0;
  calculate_boundary_conditions_call(u);
  result.push_back(u);
  alpha = dt/(dx*dx);

  for(unsigned j = 1; j < M; j++){
    u = (I-alpha*A)*u;
    calculate_boundary_conditions_call(u);
    result.push_back(u);
  }

  return result;
}


matrix FDMEulerImplicit::solve(){
    matrix result;
    BuildMatrix(); //building A and I;

    Eigen::VectorXd u = u0;
    calculate_boundary_conditions_call(u);
    result.push_back(u);
    alpha = dt/(dx*dx);

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(I + alpha*A);

    if(solver.info()!=Eigen::Success){
      std::cerr << "decomposition failed" << std::endl;
  }

    for(unsigned j = 1; j < M; j++){
      u = solver.solve(u);

      if(solver.info()!=Eigen::Success){
        std::cerr << "solving failed" << std::endl;
      }

      calculate_boundary_conditions_call(u);
      result.push_back(u);
    }

    return result;
}


matrix FDMCranckNicholson::solve(){
    matrix result;
    BuildMatrix(); //building A and I;

    Eigen::VectorXd u = u0;
    calculate_boundary_conditions_call(u);
    result.push_back(u);
    alpha = dt/(dx*dx);

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(I + alpha*A);

    if(solver.info()!=Eigen::Success){
      std::cerr << "decomposition failed" << std::endl;
  }

    for(unsigned j = 1; j < M; j++){
      u = solver.solve( (I-alpha*A) * u);

      if(solver.info()!=Eigen::Success){
        std::cerr << "solving failed" << std::endl;
      }

      calculate_boundary_conditions_call(u);
      result.push_back(u);
    }

    return result;
}
