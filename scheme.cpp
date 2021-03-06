#include "scheme.h"
#include <iomanip>
void restrict(model::Vector& vec,const model::Vector& min_corner,const model::Vector& max_corner){
  model::Vector from_left_corner =  vec - min_corner;
  double Lx = max_corner.x - min_corner.x;
  double Ly = max_corner.y - min_corner.y;
  double Lz = max_corner.z - min_corner.z;
  from_left_corner.x = Lx*(from_left_corner.x/Lx - floor(from_left_corner.x/Lx));
  from_left_corner.y = Ly*(from_left_corner.y/Ly - floor(from_left_corner.y/Ly));
  from_left_corner.z = Lz*(from_left_corner.z/Lz - floor(from_left_corner.z/Lz));
  vec = from_left_corner + min_corner;
}

Box::Box(const  std::vector<Point>& points,model::Vector min_corner_,model::Vector max_corner_):
    points(points),
    min_corner(min_corner_),
    max_corner(max_corner_),
    force(std::vector<model::Vector>(points.size())){}
VerleScheme::VerleScheme(Box& box_,int mi_,int mj_,int mk_):
box(box_),
mi(mi_),
mj(mj_),
mk(mk_),
step_number(0)
{
  model::Vector lbc = box.min_corner;
  model::Vector diag = box.max_corner-box.min_corner;
  model::Vector ex = model::Vector(diag.x,0,0);
  model::Vector ey = model::Vector(0,diag.y,0);
  model::Vector ez = model::Vector(0,0,diag.z);
  for (int i = 0;i<(2*mi+1);i++){
    left_bottom_corner.push_back(std::vector<std::vector<model::Vector>>{}); 
    for (int j = 0;j<(2*mj+1);j++){
      left_bottom_corner[i].push_back(std::vector<model::Vector>{});
      for (int k = 0;k<(2*mk+1);k++){
        model::Vector to_put;
        to_put = lbc - (mi-i)*ex - (mj-j)*ey - (mk-k)*ez;
        left_bottom_corner[i][j].push_back(to_put);
      }
    }
  }
}
model::Vector VerleScheme::GetBox(int i,int j,int k){ // returns left_bottom_corner of the box;
  return left_bottom_corner[mi-i][mj-j][mk-k];
}
/*template<class Potencial>
void VerleScheme::CalculateStep(Potencial potencial,std::vector<VerleScheme::Coord> coords){
  //std::vector<model::Vector> new_points(box.points);
  //std::vector<model::Vector> new_velocity(box.velocity);
  //std::vector<model::Vector> new_force(box.force);
  Coord main_box_coord = Coord {0,0,0};
  coords.push_back(main_box_coord);
  model::Vector diag = box.max_corner-box.min_corner;
  for (int i =0;i<box.points.size();i++){
      box.force[i] = 0;
      box.U = 0;
      for (int number = 0;number<coords.size();number++){
        model::Vector min_box_corner = left_bottom_corner[coords[number].i][coords[number].j][coords[number].k];
        if (main_box_coord==coords[number]){
            for (int j=0;j<box.points.size();j++){
              if (i!=j){
                box.force[i]+=potencial.GetForce(box.points[i],box.points[j]);
                box.U+=potencial.GetPotencial(box.points[i],box.points[j]);
              }
            }
          }
        else{
          model::Vector lbc = GetBox(coords[number].i,coords[number].j,coords[number].k);
          for (int j=0;j<box.points.size();j++){
            model::Vector new_point = box.points[j] - box.min_corner + lbc;
            box.force[i]+= potencial.GetForce(box.points[i],new_point);
            box.U+=potencial.GetPotencial(box.points[i],new_point);
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
template<class Potencial>
void VerleScheme::DoStep(double delta_t,Potencial potencial,std::vector<Coord> coords ){
  if (step_number==0)
    CalculateStep(potencial,coords);
  for (int i=0;i<box.points.size();i++){
    box.points[i]=box.points[i]+box.velocity[i]*delta_t + box.force[i]*delta_t*delta_t/(2*box.mass[i]);
  }
  std::vector<model::Vector> last_force = box.force;
  CalculateStep(potencial,coords);
  for (int i=0;i<box.velocity.size();i++){
    box.velocity[i]=box.velocity[i]+(box.force[i]+last_force[i])/(2*box.mass[i]);
  }
  step_number++;
}*/
//template<class Potencial>
/*void VerleScheme::PrintState(std::ostream& os,Potencial p){
  os<<"step"<<step_number<<"\n";
  for(int i = 0;i < box.points.size();i++){
    os<<"r"<<i+1<<"=(rx"<<i+1<<";"<<"ry"<<i+1<<";"<<"rz"<<i+1<<")=("<<std::setprecision(8)<<std::fixed<<box.points[i].x<<";"<<box.points[i].y<<";"<<box.points[i].z<<")\n";
  }
  for (int i=0;i<box.velocity.size();i++){
    os<<"v"<<i+1<<"=(vx"<<i+1<<";"<<"vy"<<i+1<<";"<<"vz"<<i+1<<")=("<<std::setprecision(8)<<std::fixed<<box.velocity[i].x<<";"<<box.velocity[i].y<<";"<<box.velocity[i].z<<")\n";
  }
  for ( int i=0;i<box.points.size()-1;i++){
    for (int j=i;j<box.points.size();j++){
      model::Vector point = box.points[j]-box.points[i];
      os<<"r"<<i+1<<j+1<<"=(rx"<<i+1<<j+1<<";"<<"ry"<<i+1<<j+1<<";"<<"rz"<<i+1<<j+1<<")=("<<std::setprecision(8)<<std::fixed<<point.x<<";"<<point.y<<";"<<point.z<<")\n";
      os<<"r"<<i+1<<j+1<<"_abs = "<<abs(point)<<"\n";

    }
  }
   os<<"U="<<box.U<<"\n";
  for ( int i=0;i<box.points.size()-1;i++){
    for (int j=i;j<box.points.size();j++){
      model::Vector f = p.GetForce(box.points[i],box.points[j]);
      os<<"F"<<i+1<<j+1<<"=(Fx"<<i+1<<j+1<<";"<<"Fy"<<i+1<<j+1<<";"<<"Fz"<<i+1<<j+1<<")=("<<std::setprecision(8)<<std::fixed<<f.x<<";"<<f.y<<";"<<f.z<<")\n";
    }
  }

}*/


