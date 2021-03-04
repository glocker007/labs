#pragma once
#include "potencial.h"
#include <cmath>
struct Point{
  model::Vector r;
  model::Vector velocity;
  double mass;
  size_t type;
};
struct Statistics{
  virtual double ChangeValue(const Point& point,const shared_ptr<Potencial::Potencial> potencial);
};
struct KineticEnergy:public Statistics{
KineticEnergy();
  double Ekin; 
  double GetValue();
  double ChangeValue(const Point& point,const shared_ptr<Potencial::Potencial> potencial) override;
};
struct KineticEnergyPerParticle:public Statistics{
KineticEnergyPerParticle();
  double GetValue();
 double Ekin; 
  size_t N;
  double ChangeValue(const Point& point,const shared_ptr<Potencial::Potencial> potencial) override; 
};
struct KineticTemperatureEnergy:public Statistics{
KineticTemperatureEnergy();
  double GetValue();
 double Etemp; 
  double ChangeValue(const Point& point,const shared_ptr<Potencial::Potencial> potencial) override;
};
struct KineticTemperatureEnergyPerParticle:public Statistics{
KineticTemperatureEnergyPerParticle();
  double GetValue();
  double Etemp;
  size_t N;
  double ChangeValue(const Point& point,const shared_ptr<Potencial::Potencial> potencial) override;
};
struct PotencialEnergy:public Statistics{
PotencialEnergy();
  double GetValue();
  double U;
  double ChangeValue(const Point& point,const shared_ptr<Potencial::Potencial> potencial) override;
};
struct PotencialEnergyPerParticle:public Statistics{
PotencialEnergyPerParticle();
  double GetValue();
  double U;
  size_t N;
  double ChangeValue(const Point& point,const shared_ptr<Potencial::Potencial> potencial) override;
};
struct FullEnergy:public Statistics{
FullEnergy();
  double GetValue();
  double E;
  double ChangeValue(const Point& point,const shared_ptr<Potencial::Potencial> potencial) override;
};
struct FullEnergyPerParticle:public Statistics{
FullEnergyPerParticle();
  double GetValue();
  double E;
  size_t N;
  double ChangeValue(const Point& point,const shared_ptr<Potencial::Potencial> potencial) override;
};
struct InternalEnergy:public Statistics{
InternalEnergy();
  double GetValue();
  double Eint;
  double ChangeValue(const Point& point,const shared_ptr<Potencial::Potencial> potencial) override;
};
struct InternalEnergyPerParticle: public Statistics{
InternalEnergyPerParticle();
  double GetValue();
  double Eint;
  size_t N;
  double ChangeValue(const Point& point,const shared_ptr<Potencial::Potencial> potencial) override;
};
struct Temperature: public Statistics{
Temperature();
  double GetValue();
  double T;
  size_t N;
  double ChangeValue(const Point& point,const shared_ptr<Potencial::Potencial> potencial) override;
};
void restrict(model::Vector& vec,const model::Vector& min_corner,const model::Vector& max_corner);
struct Box{
   /* std::vector<model::Vector> points;
    std::vector<double> mass;
    std::vector<model::Vector> velocity;*/
    std::vector<Point> points;
    std::vector<model::Vector> force;
    //std::vector<model::Vector> new_force;
    double U;
    model::Vector min_corner;
    model::Vector max_corner;
    Box(const std::vector<Point>& points_,model::Vector min_corner_,model::Vector max_corner_);
};
class VerleScheme{
  private:
    std::vector<std::vector<std::vector<model::Vector>>> left_bottom_corner;
    int mi,mj,mk;
    Box& box; 
    int step_number;
  public:
    struct Coord{
      int i;
      int j;
      int k;
      friend bool operator==(const Coord& a,const Coord& b){
        return (a.i==b.i)&&(a.j==b.j)&&(a.k==b.k);
      };
    };
    VerleScheme(Box& box,int mi,int mj,int mk);
    model::Vector GetBox(int i,int j,int k);   
    template<class Potencial_>
void CalculateStep(std::vector<std::vector<std::shared_ptr<Potencial_>>> potencial,std::vector<VerleScheme::Coord> coords){
  //std::vector<model::Vector> new_points(box.points);
  //std::vector<model::Vector> new_velocity(box.velocity);
  //std::vector<model::Vector> new_force(box.force);
  Coord main_box_coord = Coord {0,0,0};
  coords.push_back(main_box_coord);
  model::Vector diag = box.max_corner-box.min_corner;
  for (int i =0;i<box.points.size();i++){
      box.force[i] =model::Vector{0,0,0};
      box.U = 0;
      for (int number = 0;number<coords.size();number++){
        model::Vector min_box_corner = left_bottom_corner[coords[number].i][coords[number].j][coords[number].k];
        if (main_box_coord==coords[number]){
            for (int j=0;j<box.points.size();j++){
              if (i!=j){
                box.force[i]=box.force[i] + potencial[box.points[i].type][box.points[j].type]->GetForce(box.points[i].r,box.points[j].r);
                box.U= box.U + potencial[box.points[i].type][box.points[j].type]->GetPotencial(box.points[i].r,box.points[j].r);
              }
            }
          }
        else{
          model::Vector lbc = GetBox(coords[number].i,coords[number].j,coords[number].k);
          for (int j=0;j<box.points.size();j++){
            model::Vector new_point = box.points[j].r - box.min_corner + lbc;
            box.force[i]= box.force[i] + potencial[box.points[i].type][box.points[i].type]->GetForce(box.points[i].r,new_point);
            box.U+=potencial[box.points[i].type][box.points[i].type]->GetPotencial(box.points[i].r,new_point);
          }
        }
      }
     // new_points[i]+=box.velocity[i]*delta_t+box.force[i]*delta_t*delta_t/(2*box.mass[i]);
      //new_velocity[i]+=(box.force[i]+new_force[i])*delta_t/box.mass[i];
  }
  //box.points = new_points;
  //box.velocity = new_velocity;
  //box.force = new_force;
}
   template<class Potencial_>
  void DoStep(double delta_t,std::vector<std::vector<std::shared_ptr<Potencial_>>> potencial,std::vector<Coord> coords ){
    if (step_number==0)
      CalculateStep(potencial,coords);
    for (int i=0;i<box.points.size();i++){
      box.points[i].r=box.points[i].r+box.points[i].velocity*delta_t + box.force[i]*delta_t*delta_t/(2*box.points[i].mass);
      restrict(box.points[i].r,box.min_corner,box.max_corner);
    }
    std::vector<model::Vector> last_force = box.force;
    CalculateStep(potencial,coords);
    for (int i=0;i<box.points.size();i++){
      box.points[i].velocity=box.points[i].velocity+(box.force[i]+last_force[i])*delta_t/(2*box.points[i].mass);
    }
    step_number++;
  }
   
   template<class Potencial_>
  void PrintState(std::ostream& os,std::vector<std::vector<std::shared_ptr<Potencial_>>> p){
  os<<"Step = "<<step_number<<"\n";
  for(int i = 0;i < box.points.size();i++){
    os<<"r"<<i+1<<"=(rx"<<i+1<<";"<<"ry"<<i+1<<";"<<"rz"<<i+1<<")=("<<std::setprecision(8)<<std::fixed<<box.points[i].r.x<<";"<<box.points[i].r.y<<";"<<box.points[i].r.z<<")\n";
  }
  for (int i=0;i<box.points.size();i++){
    os<<"v"<<i+1<<"=(vx"<<i+1<<";"<<"vy"<<i+1<<";"<<"vz"<<i+1<<")=("<<std::setprecision(8)<<std::fixed<<box.points[i].velocity.x<<";"<<box.points[i].velocity.y<<";"<<box.points[i].velocity.z<<")\n";
  }
  for ( int i=0;i<box.points.size()-1;i++){
    for (int j=i+1;j<box.points.size();j++){
      model::Vector point = box.points[j].r-box.points[i].r;
      os<<"r"<<i+1<<j+1<<"=(rx"<<i+1<<j+1<<";"<<"ry"<<i+1<<j+1<<";"<<"rz"<<i+1<<j+1<<")=("<<std::setprecision(8)<<std::fixed<<point.x<<";"<<point.y<<";"<<point.z<<")\n";
      os<<"r"<<i+1<<j+1<<"_abs = "<<abs(point)<<"\n";

    }
  }
   os<<"U="<<box.U<<"\n";
  for ( int i=0;i<box.points.size()-1;i++){
    for (int j=i+1;j<box.points.size();j++){
      model::Vector f = p[box.points[i].type][box.points[j].type]->GetForce(box.points[i].r,box.points[j].r);
      os<<"F"<<i+1<<j+1<<"=(Fx"<<i+1<<j+1<<";"<<"Fy"<<i+1<<j+1<<";"<<"Fz"<<i+1<<j+1<<")=("<<std::setprecision(24)<<std::fixed<<f.x<<";"<<f.y<<";"<<f.z<<")\n";
    }
  }

    }
};
