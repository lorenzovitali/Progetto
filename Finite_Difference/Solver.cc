#include "Solver.hh"

void Solver::set_initial_conditions() {
  double cur_spot = Nminus*dx;
  u0.resize(N+1);
  x_values.reserve(N+1);

  for (unsigned j = 0; j <= N; j++) {
    u0[j] = pde->init_cond(cur_spot);
    x_values.push_back(cur_spot);
    cur_spot += dx;
  }

  tau_values.reserve(M+1);
  double time_step = 0.0;
  for (unsigned i = 0; i <= M; ++i){
    tau_values.push_back(time_step);
    //std::cout << "tau: " << tau_values[i] << std::endl;
    time_step += dt;
  }
}

void Solver::calculate_boundary_conditions(std::vector<double> & u_new, double tau) {
  u_new[0] = pde->boundary_left(Nminus*dx, tau);
  u_new[N] = pde->boundary_right(Nplus*dx, tau);
}

void Solver::change_var(Matrice& matrice){
    for(unsigned i = 0; i <= M; ++i){
      for(unsigned j = 0; j <= N; ++j){
        matrice[i][j] = pde->get_option()->get_E()*matrice[i][j] * exp(-0.5 * (pde->get_option()->get_k() - 1) * x_values[j] - 0.25*(pde->get_option()->get_k() + 1)*(pde->get_option()->get_k() + 1) *tau_values[i]);
      }
    }
  }

double Solver::find_opt(){
  std::vector<double> omega;
  unsigned size = 100;
  double current = 1.0;
  double step = 1/static_cast<double>(size);

  for(unsigned i = 0; i < size; ++i){
    omega.push_back(current);
    current += step;
  }

  std::vector<unsigned> iter_vec;
  for(unsigned i = 0; i < omega.size(); ++i){
    solve(omega[i]);
    iter_vec.push_back(iter);
  }

  unsigned index = 0;
  for(unsigned i = 1; i < iter_vec.size(); ++i){
    if(iter_vec[i] < iter_vec[i-1]) index = i;
  }

  std::cout << "omega optimal: " << omega[index] << std::endl;
  std::cout << "iter: " << iter_vec[index] << std::endl;

  return omega[index];
}


Matrice SOR_EI::solve(double w){
  Matrice risultato;
  set_initial_conditions();
  std::vector<double> u = u0;

  calculate_boundary_conditions(u,tau_values[0]);
  //std::cout << pde->boundary_right(Nplus*dx, tau_values[0]) << std::endl;
  risultato.push_back(u);
  double err;
  alpha = dt/(dx*dx);
  double y;
  std::vector<double> b;

  iter = 0;
  for(unsigned i = 1; i <= M; ++i){

    b = risultato[i-1];

    b[1] += alpha*pde->boundary_left(Nminus*dx, tau_values[i]);
    b[N-1] += alpha*pde->boundary_right(Nplus*dx, tau_values[i]);

    //std::cout << pde->boundary_right(Nplus*dx, tau_values[i]) << std::endl;

    err = 10000;

    u = risultato[i-1];

    while( err > toll ){

      err = 0.0;

      for(unsigned j = 1; j < N ; ++j){

        //std::cout << "x: " << x_values[j] << std::endl;
        y = ( b[j] + alpha*(u[j-1] + u[j+1]) )/(1+2*alpha);
        y = u[j] + w*(y - u[j]);
        err += (u[j] - y )*(u[j] - y);
        u[j] = y;

      }
      iter++;
      //std::cout << "err: " << err << std::endl;
    }
    calculate_boundary_conditions(u,tau_values[i]);
    risultato.push_back(u);
  }
  change_var(risultato);
  return risultato;
}


Matrice PSOR_CN::solve(double w){
  Matrice risultato;
  set_initial_conditions();
  std::vector<double> u = u0;

  calculate_boundary_conditions(u,tau_values[0]);
  //std::cout << pde->boundary_right(Nplus*dx, tau_values[0]) << std::endl;
  risultato.push_back(u);
  double err;
  alpha = dt/(dx*dx);
  double y;
  double b;
  double g;

  iter = 0;
  for(unsigned i = 1; i <= M; ++i){

    err = 10000;

    u = risultato[i-1];

    while( err > toll ){

      err = 0.0;

      for(unsigned j = 1; j < N ; ++j){

        if(j==1){
          b = (1 - alpha)*risultato[i-1][j] + alpha/2*(risultato[i-1][j+1] + risultato[i-1][j-1]) + alpha/2*pde->boundary_left(Nminus*dx, tau_values[i]);
        }
        else if(j==N-1)
            b = (1 - alpha)*risultato[i-1][j] + alpha/2*(risultato[i-1][j+1] + risultato[i-1][j-1]) + alpha/2*pde->boundary_right(Nplus*dx, tau_values[i]);
        else{
            b = (1 - alpha)*risultato[i-1][j] + alpha/2*(risultato[i-1][j+1] + risultato[i-1][j-1]);
        }
        //std::cout << "x: " << x_values[j] << std::endl;

        y = ( b + alpha/2*(u[j-1] + u[j+1]) )/(1+alpha);
        g = exp(0.25 * (pde->get_option()->get_k()+1)*(pde->get_option()->get_k()+1) * tau_values[i]) * pde->init_cond(x_values[j]);
        y = std::max( u[j] + w*(y - u[j]), g);

        err += (u[j] - y )*(u[j] - y);
        u[j] = y;

      }
      iter++;
      //std::cout << "err: " << err << std::endl;
    }
    calculate_boundary_conditions(u,tau_values[i]);
    risultato.push_back(u);
  }
  change_var(risultato);
  return risultato;
}
