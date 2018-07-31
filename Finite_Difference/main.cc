#include "Option.hh"
#include "PDE.hh"
#include "FDM.hh"
#include <iostream>

int main(){
  Option* call = new EuropeanCall(0.5, 0.05, 0.2, 1.00);
  BlackScholesPDE* pde = new BlackScholesPDE(call);

  FDMEulerExplicit fdm_EE (20,10, pde);



  /*for(unsigned i=0;i<10 ;++i){
    for(unsigned j=0;j<20;++j){
      std::cout << solution[i][j] << "\t";
    }
    std::cout << std::endl;
  }*/



  return 0;
}
