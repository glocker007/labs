#include "potencial.h"
#include <cmath>
namespace Potencial{
  LenardJones::LenardJones(int m_,int n_,double sigma_,double epsilon_):
    m(m_),
    n(n_),
    sigma(sigma_),
    epsilon(epsilon_)
  {
  }
  model::Vector LenardJones::GetForce(const model::Vector& v1,const model::Vector& v2){
    double distance = abs(v1-v2);
    return 4*epsilon*(m * (pow(sigma,m)/pow(distance,m+2)) * (v1-v2) - n * (pow(sigma,n)/pow(distance,n+2)) * (v1-v2));
  }
  double LenardJones::GetPotencial(const model::Vector& v1,const model::Vector& v2){
    double distance = abs(v1-v2);
    return 4*epsilon*(pow(sigma,12)/pow(distance,12) - pow(sigma,6)/pow(distance,6));
  }
}


