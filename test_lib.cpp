#include "scheme.h"
#include <memory>
using namespace std;
int main(){
  vector<Point> points;
  points.push_back(Point{model::Vector(0.25,0.75,0.5),model::Vector(1,1,0),66.335,0});
  points.push_back(Point{model::Vector(0.75,0.75,0.5),model::Vector(-1,1,0),66.335,0});
  model::Vector min_corner = model::Vector{0,0,0};
  model::Vector max_corner = model::Vector(5,5,5);
  vector<model::Vector> force;
  force.resize(points.size());
  Box box(points,min_corner,max_corner);
  Potencial::LenardJones lj(12,6,0.3418,1.712);
  vector<vector<shared_ptr<Potencial::Potencial>>> vec; 
  vec.push_back({make_shared<Potencial::LenardJones>(lj)});
  VerleScheme  scheme(box,0,0,0);
  model::Vector v10(0.25200020,0.75200000,0.50000000);
  model::Vector v20(0.74799980,0.75200000,0.50000000);
  model::Vector v11(0.252000201234882537093540,0.752000000000000001776357,0.500000000000000000000000);
  model::Vector v21(0.747999798765117462906460,0.752000000000000001776357,0.500000000000000000000000);
  cout<<lj.GetForce(v10,v20)<<endl;
  cout<< lj.GetForce(v11,v21)<<endl;  
  scheme.CalculateStep(vec,{});
  scheme.PrintState(cout,vec);
  cout<<"\n";
  scheme.DoStep(0.002,vec,{});
  scheme.PrintState(cout,vec);

  return 0;
}
