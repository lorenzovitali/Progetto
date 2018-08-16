#ifndef PDE_HH
#define PDE_HH

#include <memory>
#include "Option.hh"

class BlackScholesPDE {
  protected:
    std::shared_ptr<Option> option;

  public:
    BlackScholesPDE(std::shared_ptr<Option> _option): option(_option) {}

    std::shared_ptr<Option> get_option()const{return option;}

    //boundaries
    virtual double boundary_left(double x, double tau) const = 0;

    virtual double boundary_right(double x, double tau) const = 0;

    double init_cond(double x);
};

class BlackScholesPDECall:public BlackScholesPDE{
  public:
    BlackScholesPDECall(std::shared_ptr<Option> _option): BlackScholesPDE(_option){}

    double boundary_left(double x, double tau) const;

    double boundary_right(double x, double tau) const;
};

class BlackScholesPDEPut:public BlackScholesPDE{
  public:
    BlackScholesPDEPut(std::shared_ptr<Option> _option): BlackScholesPDE(_option){}

    double boundary_left(double x, double tau) const;

    double boundary_right(double x, double tau) const;
};
#endif
