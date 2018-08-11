#include "SOR.hh"

void SOR::set_initial_conditions() {
  double cur_spot = Nminus*dx;
  u0.resize(N);
  x_values.reserve(N);

  for (unsigned j = 0; j < N; j++) {
    u0[j] = pde->init_cond(cur_spot);
    x_values.push_back(cur_spot);
    cur_spot += dx;
  }
}

void SOR::calculate_boundary_conditions(std::vector<double> & u_new, double tau) {
  u_new[0] = pde->boundary_left(Nminus*dx, tau);
  u_new[N-1] = pde->boundary_right(Nplus*dx, tau);
}

void SOR::change_var(Matrice& matrice){
    for(unsigned i = 0; i < M; ++i){
      for(unsigned j = 0; j < N; ++j){
        matrice[i][j] = pde->get_option()->get_E()*matrice[i][j] * exp(-0.5 * (pde->get_option()->get_k() - 1) * x_values[j] - 0.25*(pde->get_option()->get_k() + 1)*(pde->get_option()->get_k() + 1) *tau_values[i]);
      }
    }
  }

Matrice SOR::solve(double w){
  Matrice risultato;
  set_initial_conditions();
  std::vector<double> u = u0;

  tau_values.reserve(M);
  double time_step = 0.0;
  tau_values.push_back(time_step);

  calculate_boundary_conditions(u,time_step);
  risultato.push_back(u);
  double err;
  alpha = dt/(dx*dx);
  double y;
  std::vector<double> b;

  iter = 0;
  for(unsigned i = 1; i < M; ++i){

    time_step += dt;
    tau_values.push_back(time_step);
    b = risultato[i-1];
    err = 10000;
    //for(unsigned index = 0; index < N; ++index) u[index] = 0;
    u = risultato[i-1];

    while( err > toll ){

      err = 0.0;

      for(unsigned j = 1; j < N - 1 ; ++j){

        y = ( b[j] + alpha*(u[j-1] + u[j+1]) )/(1+2*alpha);
        y = u[j] + w*(y - u[j]);
        err += (u[j] -y )*(u[j] - y);
        u[j] = y;

      }
      iter++;
      //std::cout << "err: " << err << std::endl;
    }
    calculate_boundary_conditions(u,time_step);
    risultato.push_back(u);
  }
  std::cout << "rows = " <<  risultato.size() << std::endl;
  std::cout << "columns = " << risultato[1].size() << std::endl;
  change_var(risultato);
  return risultato;
}
