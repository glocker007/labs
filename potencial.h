#pragma once
#include "model_vector.h"
namespace Potencial{
  struct Potencial{
    int m = 0 ;
    virtual model::Vector GetForce(const model::Vector& v1,const model::Vector& v2)  =  0;
    virtual double GetPotencial(const model::Vector& v1, const model::Vector& v2) = 0;
  };
  struct LenardJones:public Potencial{
    int m;
    int n;
    double sigma;
    double epsilon;
    LenardJones() = default;
    LenardJones(int m,int n,double sigma,double epsilon); 
    model::Vector GetForce(const model::Vector& v1,const model::Vector& v2) override;
    double GetPotencial(const model::Vector& v1,const model::Vector& v2) override;
  };
}

