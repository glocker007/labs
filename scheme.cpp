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



