#include "model_vector.h"
#include <string>
namespace model{
  Vector::Vector(){
    x=0;
    y=0;
    z = 0;
  }
Vector::Vector(double x_,double y_,double z_):x(x_),y(y_),z(z_){}


  Vector& Vector::operator=(const Vector& vec){
    x = vec.x;
    y = vec.y;
    z = vec.z;
    return *this;
  }
  //Vector& operator=(const MaterialPoint& mp);
  
  MaterialPoint::MaterialPoint(double x_,double y_,double z_,double mass_):
    Vector(x_,y_,z_),mass(mass_)
  {
  }
}

double s(const model::Vector& v1,const model::Vector& v2){
  return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;
}

  
model::Vector operator*(double a,const model::Vector& v){
  return model::Vector{a*v.x,a*v.y,a*v.z};
}
model::Vector operator*(const model::Vector& v,double a){
  return model::Vector{a*v.x,a*v.y,a*v.z};
}
model::Vector operator/(const model::Vector& v,double a){
  return model::Vector{v.x/a,v.y/a,v.z/a};
}

model::Vector operator-(const model::Vector& v1,const model::Vector& v2){
  return model::Vector{v1.x-v2.x,v1.y-v2.y,v1.z-v2.z};
}

model::Vector operator+(const model::Vector& v1,const model::Vector& v2){
  return model::Vector{v1.x+v2.x,v1.y+v2.y,v1.z+v2.z};
}

double abs(model::Vector v1,model::Vector v2){
  return pow(s((v1-v2),(v1-v2)),0.5);
}

double abs(model::Vector v1){
  return pow(s(v1,v1),0.5);
}
std::ostream& operator<<(std::ostream& os,const model::Vector& v1) {
    os<<std::string("(")<<std::setprecision(9)<<v1.x<<std::string("; ")<<v1.y<<std::string("; ")<<v1.z<<std::string(")");
    return os;
  }

