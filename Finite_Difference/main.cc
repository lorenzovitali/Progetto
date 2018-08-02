#include "Option.hh"
#include "PDE.hh"
#include "FDM.hh"
#include <iostream>

int main(){
  double E = 0.5;
  double r = 0.05;
  double T = 1.00;
  double sigma = 0.2;

  Option* call = new EuropeanCall(E, r, T, sigma);
  BlackScholesPDE* pde = new BlackScholesPDE(call);

  FDMEulerExplicit fdm_EE (20,10, pde);
  matrix solution;
  solution = fdm_EE.solve() ;


  for(unsigned i=0;i<10 ;++i){
    for(unsigned j=0;j<20;++j){
      std::cout << solution[i][j] << "\t";
    }
    std::cout << std::endl;
  }



  return 0;
}
