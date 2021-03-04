#pragma once
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <memory>
namespace model{
struct Vector{
  double x;
  double y;
  double z;  

  Vector();
  
  Vector(double x_,double y_,double z_);
  

  
  Vector& operator=(const Vector& vec);
  //Vector& operator=(const MaterialPoint& mp);
};
struct MaterialPoint: Vector{
  double mass;
  
  MaterialPoint(double x_,double y_,double z_,double mass_);
};
}
double s(const model::Vector& v1,const model::Vector& v2);

  
model::Vector operator*(double a,const model::Vector& v);
model::Vector operator*(const model::Vector& v,double a);
model::Vector operator/(const model::Vector& v,double a);

model::Vector operator-(const model::Vector& v1,const model::Vector& v2);

model::Vector operator+(const model::Vector& v1,const model::Vector& v2);

double abs(model::Vector v1,model::Vector v2);

double abs(model::Vector v1);

std::ostream& operator<<(std::ostream& os,const model::Vector& v1);
