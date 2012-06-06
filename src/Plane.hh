#ifndef PLANE_HH
#define PLANE_HH
#include <vector>
class Plane
{
  public:
    Plane(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax)
    {
       min.resize(3);
       max.resize(3);
       min[0] = xmin;
       min[1] = ymin;
       min[2] = zmin;
       max[0] = xmax;
       max[1] = ymax;
       max[2] = zmax;
    };
    Plane()
    {
       min.resize(3);
       max.resize(3);
    };
    ~Plane() {};
    std::vector<int> min;
    std::vector<int> max;
    
private:

};
#endif
