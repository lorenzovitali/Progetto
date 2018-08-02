#include "Option.hh"
#include "PDE.hh"
#include "FDM.hh"
#include "GetPot"
#include <iostream>

int main(int argc, char* argv[]){
  
  GetPot file("input.txt");
  double E = file("E", 0.0);
  double r = file("r", 0.05);
  double T = file("T", 1.00);
  double sigma = file("sigma", 0.2);
  unsigned n = file("n", 1000); //space intervals
  unsigned M = file("M", 10); //time intervals

  Option* call = new EuropeanCall(E, r, T, sigma);
  BlackScholesPDE* pde = new BlackScholesPDE(call);

  FDMEulerExplicit fdm_EE (n, M, pde);
  matrix solution;
  solution = fdm_EE.solve() ;

  for(unsigned i=0;i<M ;++i){
    for(unsigned j=0;j<n;++j){
      std::cout << solution[i][j] << "\t";
    }
    std::cout << "\n" << std::endl;
  }

  std::cout << "alpha: " << fdm_EE.get_alpha() << std::endl;
  std::cout << "dx: " << fdm_EE.get_dx() << std::endl;
  std::cout << "dt: "<< fdm_EE.get_dt() << std::endl;

  std::cout << "rows: "<< solution.size() << std::endl;
  std::cout << "columns: " << solution[0].size() << std::endl;


  return 0;
}
