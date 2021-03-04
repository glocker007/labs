#include <iostream>
#include "scheme.h"
using namespace std;
int main(){
  model::Vector v1(-7.2,-3.4,5.8);
  model::Vector min_corner(1,1,1);
  model::Vector max_corner(3,4,5);
  restrict(v1,min_corner,max_corner);
  cout<<v1<<endl;
  return 0;
}
