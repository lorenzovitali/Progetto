#ifndef SOR_HH
#define SOR_HH
#include "PDE.hh"
#include <vector>

typedef std::vector<std::vector<double>> Matrice;

class SOR{
  private:
     std::shared_ptr<BlackScholesPDE> pde;
     std::vector<double> u0;
     unsigned iter;

     double alpha;
     double toll;

     int long N;
     int long Nminus;
     unsigned Nplus;
     double dx;
     std::vector<double> x_values;

     std::vector<double> tau_values;
     unsigned long M;
     double dt;

  public:
    SOR(int long _N, double _dx, unsigned long _M, double _toll, std::shared_ptr<BlackScholesPDE> _pde):
      N(_N), dx(_dx), M(_M), toll(_toll), pde(_pde){
        Nminus = -N/2;
        Nplus = N/2;
        dt = (0.5* pde->get_option()->get_sigma() * pde->get_option()->get_sigma() * pde->get_option()->get_T())/M;
      }

    Matrice solve(double w);

    void set_initial_conditions();

    void calculate_boundary_conditions(std::vector<double>& , double);

    std::vector<double> get_x() const {return x_values; }
    std::vector<double> get_tau() const {return tau_values; }
    unsigned get_iter()const{return iter;}

    void change_var(Matrice&);

    //double find_opt() const;
};

#endif
