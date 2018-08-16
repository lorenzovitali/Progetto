#include "Option.hh"
#include "PDE.hh"
#include "FDM.hh"
#include "GetPot"
#include "SOR.hh"
#include <iostream>

int main(int argc, char* argv[]){

  GetPot file("input.txt");
  double E = file("E", 0.0);
  double r = file("r", 0.05);
  double T = file("T", 1.00);
  double sigma = file("sigma", 0.2);
  int long N = file("N", 100); //space intervals
  double toll = file("toll", 0.001);
  double dx = file("dx", 0.05);
  unsigned long M = file("M", 10); //time intervals

  std::shared_ptr<Option> call = std::make_shared<EuropeanCall> (EuropeanCall(E, r, T, sigma));
  std::shared_ptr<BlackScholesPDE> pde = std::make_shared<BlackScholesPDECall> (BlackScholesPDECall(call));

  FDMEulerExplicit fdm (N, dx, M, pde);
  matrix solution;
  solution = fdm.solve() ;

  /*SOR solver(N, dx, M, toll, pde);
  Matrice solution = solver.solve(1.5);*/

  for(unsigned i = 0; i < M ;++i){
    for(unsigned j = 0; j < N; ++j){
      std::cout << solution[i][j] << "\t";
    }
    std::cout << "\n" << std::endl;
  }

  //std::cout << "iter: " << solver.get_iter() << std::endl;

  /*std::cout << "alpha: " << fdm.get_alpha() << std::endl;
  std::cout << "dx: " << fdm.get_dx() << std::endl;
  std::cout << "dt: "<< fdm.get_dt() << std::endl;
  std::cout << "k: " << call->get_k() << std::endl;

  std::cout << "rows: "<< solution.size() << std::endl;
  std::cout << "columns: " << solution[0].size() << std::endl;*/


  return 0;
}
