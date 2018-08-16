#include "Option.hh"
#include "PDE.hh"
#include "FDM.hh"
#include "GetPot"
#include "Solver.hh"
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



  FDMEulerExplicit fdm_EE (N, dx, M, pde);
  matrix result;
  result = fdm_EE.solve() ;

  std::cout << "alpha: " << fdm_EE.get_alpha() << std::endl;
  std::cout << "dx: " << fdm_EE.get_dx() << std::endl;
  std::cout << "dt: "<< fdm_EE.get_dt() << std::endl;
  std::cout << "k: " << call->get_k() << std::endl;

  std::cout << " \n\n\n EULERO ESPLICITO \n" << std::endl;
  for(unsigned i = 0; i <= M ;++i){
    for(unsigned j = 0; j <= N; ++j){
      std::cout << result[i][j] << "\t";
    }
    std::cout << "\n" << std::endl;
  }

  std::cout << "rows: "<< result.size() << std::endl;
  std::cout << "columns: " << result[0].size() << std::endl;

  FDMEulerImplicit fdm_EI(N, dx, M, pde);
  matrix solution = fdm_EI.solve();

  std::cout << "\n \n \nEULERO IMPLICITO \n" << std::endl;
  for(unsigned i = 0; i <= M ;++i){
    for(unsigned j = 0; j <= N; ++j){
      std::cout << solution[i][j] << "\t";
    }
    std::cout << "\n" << std::endl;
  }

  FDMCranckNicholson fdm_CN(N, dx, M, pde);
  matrix soluzione = fdm_CN.solve();

  std::cout << "\n \n \nCRANCK NICHOLSON \n" << std::endl;
  for(unsigned i = 0; i <= M ;++i){
    for(unsigned j = 0; j <= N; ++j){
      std::cout << soluzione[i][j] << "\t";
    }
    std::cout << "\n" << std::endl;
  }


  /*PSOR_CN solver(N, dx, M, toll, pde);
  //double omega = solver.find_opt();
  Matrice risultato = solver.solve(1.04);

  for(unsigned i = 0; i <= M ;++i){
    for(unsigned j = 0; j <= N; ++j){
      std::cout << risultato[i][j] << "\t";
    }
    std::cout << "\n" << std::endl;
  }*/


  return 0;
}
