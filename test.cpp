#include <iostream>
#include "potencial.h"

using namespace std;
int main(){
  double a;
  Potencial::LenardJones lj(12,6,0.3418,1.712);
  model::Vector v1(0.25,0.75,0.5);
  model::Vector v2(0.75,0.75,0.5);
  cout<<lj.GetPotencial(v1,v2)<<endl;
  cout<<lj.GetForce(v1,v2)<<endl;
  cout<<lj.GetForce(v2,v1)<<endl;
  model::Vector vec = {-0.7,2.81,9.2};
  //cout<<bounds(vec,Vector(-0.5,0.2,0.5)) <<endl;
  cout<<vec*2-vec+1.2*vec/5<<endl;
  return 0;
}
