#ifndef PROC_BOX_HH
#define PROC_BOX_HH

#include <set>
#include <map>
#include <vector>
#include "Plane.hh"

class ProcBox
{
  public:
    ProcBox(int peid);
    ~ProcBox();
    int volume();
    int myRank() { return peid_; }
    void add3DCoord(int x, int y, int z);
    int nAdded() { return addCnt_; }
    int trialMove(int dim, int val, ProcBox& nbrbox);
    void planesToCoords(std::vector<int>& coords);
    void addPlane(int dim, Plane* plptr);
    void removePlane(int dim, int val);
    void deletePlane(int dim, int val);
    void updateBoundaries();
    Plane* getPlane(int dim, int val) { return planes_[dim][val]; }
    int minAxisPoint(int dim);
    int maxAxisPoint(int dim);
    int secondMinAxisPoint(int dim);
    int secondMaxAxisPoint(int dim);
    int minPoint(int dim) { return min_[dim]; }
    int maxPoint(int dim) { return max_[dim]; }
    bool containsPoint(int x, int y, int z);
    int trialVolume(bool removePlane, int dim, Plane* pptr);
    
private:
    int peid_;
    int addCnt_;
    std::vector<std::map<int,Plane*> > planes_;
    std::vector<std::set<int> > axisPoints_;
    std::vector<int> min_;
    std::vector<int> max_;
};


#endif
