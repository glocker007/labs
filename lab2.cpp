#include <iostream>
#include "scheme.h"
#include <fstream>
using namespace std;
int main(){
  fstream file("Ganiev_lab2.txt");
  fstream file2("for_graph.txt");
  fstream file3("for_graph2.txt");
  vector<Point> points;
  points.push_back(Point{model::Vector(0.25,0.75,0.5),model::Vector(1,1,0),66.335,0});
  points.push_back(Point{model::Vector(0.75,0.75,0.5),model::Vector(-1,-1,0),66.335,0});
  points.push_back(Point{model::Vector(1.35,1.75,0.5),model::Vector(1,2,0),66.335,0});
  points.push_back(Point{model::Vector(1.68,1.75,0.5),model::Vector(-1,-2,0),66.335,0});
  Potencial::LenardJones lj(12,6,0.3418,1.712);//k = 1.38065
  vector<vector<shared_ptr<Potencial::Potencial>>> vec;
  vec.push_back({make_shared<Potencial::LenardJones>(lj)});
  Box b(points,model::Vector{0,0,0},model::Vector{4,4,4});
  VerleScheme vs(b,0,0,0);
  Stat avg;
  for (int i=0;i<2000;i++){
      vs.PrintCoords(file2);
      Stat s = vs.DoStep(0.002,vec,{});
      for (int j = 0;j<6;j++){
          Statistics a = static_cast<Statistics>(j);
          if (j!=5) file3<<fixed<<setprecision(8)<<s[a]/4.<<" ";
          else file3<<fixed<<setprecision(8)<<s[a];
      }
      file3<<"\n";
      vs.PrintStats(file,i);
      avg=avg + (1./(2000)) * s; 
  }
  avg = 0.25*avg;
  avg[Statistics::T]*=4;
  cout<<avg<<endl; 

  return 0;
}

