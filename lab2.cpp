#include <iostream>
#include "scheme.h"
using namespace std;
int main(){
  vector<Points> points;
  points.push_back(Point{model::Vector(0.25,0.75,0.5),model::Vector(1,1,0),66.335,0});
  points.push_back(Point{model::Vector(0.75,0.75,0.5),model::Vector(-1,-1,0),66.335,0});
  points.push_back(Point{model::Vector(1.35,1.75,0.5),model::Vector(1,2,0),66.335,0});
  points.push_back(Point{model::Vector(1.68,1.75,0.5),model::Vector(-1,-1,0),66.335,0});
  Potencial::LenardJones lj(12,6,0.3418,1.712);
  return 0;
}

