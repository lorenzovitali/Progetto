#include <iostream>
#include "FDM.hh"

FDMBase::FDMBase(int _N, double _dx, unsigned long _M, std::shared_ptr<BlackScholesPDE> _pde) :N(_N), dx(_dx), M(_M), pde(_pde){

                                      calculate_step_sizes();
                                      set_initial_conditions();

                                    }

void FDMBase::BuildMatrix(){
  //bulding matrix A and I;
  A.resize(N+1,N+1);
  I.resize(N+1,N+1);

  using T = Eigen::Triplet<double>;
  std::vector<T> coeff;

  coeff.emplace_back(0, 0, 2);
  coeff.emplace_back(0, 1, -1);
  coeff.emplace_back(N , N, -1);
  coeff.emplace_back(N ,N-1, 2);

  for (unsigned i = 1; i <= N-1; i++){
    coeff.emplace_back(i,i,2);
    coeff.emplace_back(i,i-1,-1);
    coeff.emplace_back(i,i+1,-1);
  }

  A.setFromTriplets(coeff.begin(), coeff.end());
  coeff.clear();

  for (unsigned i = 0; i <= N; i++){
    coeff.emplace_back(i,i,1);
  }

  I.setFromTriplets(coeff.begin(), coeff.end());
}


void FDMBase::calculate_step_sizes() {

  Nplus = N/2;
  Nminus = -N/2;

  dt = (0.5* pde->get_option()->get_sigma() * pde->get_option()->get_sigma() * pde->get_option()->get_T())/M;

}


void FDMBase::set_initial_conditions() {
  // Spatial settings
  double cur_spot = Nminus*dx;
  x_values.reserve(N+1);
  u0.resize(N+1);

  for (unsigned long j = 0; j <= N; j++) {

    u0[j] = pde->init_cond(cur_spot);
    x_values.push_back(cur_spot);
    //std::cout << "x: " << x_values[j] << std::endl;
    cur_spot += dx;

  }

  double time_step = 0.0;
  tau_values.reserve(M+1);
  for (unsigned i = 0; i <= M; ++i){
    tau_values.push_back(time_step);
    //std::cout << "tau: " << tau_values[i] << std::endl;
    time_step += dt;
  }
}


void FDMBase::calculate_boundary_conditions(Eigen::VectorXd & u_new, double tau) {
  u_new[0] = pde->boundary_left(Nminus*dx, tau);
  u_new[N] = pde->boundary_right(Nplus*dx, tau);
}

void FDMBase::change_var(matrix& matrice){
    for(unsigned i = 0; i <= M; ++i){
      for(unsigned j = 0; j <= N; ++j){
        matrice[i][j] = pde->get_option()->get_E()*matrice[i][j] * exp(-0.5 * (pde->get_option()->get_k() - 1) * x_values[j] - 0.25*(pde->get_option()->get_k() + 1)*(pde->get_option()->get_k() + 1) *tau_values[i]);
      }
    }
  }

matrix FDMEulerExplicit::solve(){
  matrix result;
  BuildMatrix(); //building A and I;
  Eigen::VectorXd u = u0;

  calculate_boundary_conditions(u, tau_values[0]);
  //std::cout << pde->boundary_right(Nplus*dx, tau_values[0]) << std::endl;
  result.push_back(u);
  alpha = dt/(dx*dx);

  for(unsigned i = 1; i <= M; i++){
    u = (I-alpha*A)*u;
    calculate_boundary_conditions(u, tau_values[i]);
    result.push_back(u);
    //std::cout << pde->boundary_right(Nplus*dx, tau_values[i]) << std::endl;
  }

  change_var(result);
  return result;
}


matrix FDMEulerImplicit::solve(){
    matrix result;
    BuildMatrix(); //building A and I;
    Eigen::VectorXd u = u0;

    calculate_boundary_conditions(u, tau_values[0]);
    result.push_back(u);
    alpha = dt/(dx*dx);

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(I + alpha*A);

    if(solver.info()!=Eigen::Success){
      std::cerr << "decomposition failed" << std::endl;
  }

    for(unsigned i = 1; i <= M; i++){
      u = solver.solve(u);

      if(solver.info()!=Eigen::Success){
        std::cerr << "solving failed" << std::endl;
      }

      calculate_boundary_conditions(u, tau_values[i]);
      result.push_back(u);
    }

    change_var(result);
    return result;
}


matrix FDMCranckNicholson::solve(){
    matrix result;
    BuildMatrix(); //building A and I;
    Eigen::VectorXd u = u0;

    calculate_boundary_conditions(u, tau_values[0]);
    result.push_back(u);
    alpha = dt/(dx*dx);

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(I + alpha*A);

    if(solver.info()!=Eigen::Success){
      std::cerr << "decomposition failed" << std::endl;
  }

    for(unsigned i = 1; i <= M; i++){
      u = solver.solve( (I-alpha*A) * u);

      if(solver.info()!=Eigen::Success){
        std::cerr << "solving failed" << std::endl;
      }

      calculate_boundary_conditions(u, tau_values[i]);
      result.push_back(u);
    }

    change_var(result);
    return result;
}
