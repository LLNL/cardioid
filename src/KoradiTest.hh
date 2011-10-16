#ifndef KORADI_TEST_HH
#define KORADI_TEST_HH

#include <vector>
#include "AnatomyCell.hh"
#include "three_algebra.h"
#include "IndexToVector.hh"

class Simulate;

class KoradiTest
{
 public:
   
   KoradiTest(Simulate& sim, int nCentersPerTask,
	      double alpha, int voronoiSteps);
   
   void balanceStep();
   void recondition(int nSteps);
   
 private:

   void distributeCellsEvenly();
   void pickInitialCenters();
   void assignCells();
   void initialAssignment();
   void sendCellsToDestinations();
   void moveCenters();
   void computeRadii();
   void printStatistics();
   void updateBias();
   
   void calculateCellDestinations();
   void exchangeCells();
   void bruteForceDistanceCheck();
   void findNbrDomains();
   void computeLoad(std::vector<double>& load);


   double _alpha;
   
   int _myRank;
   int _nTasks;

   int _nCentersPerTask;
   int _localOffset;

   IndexToVector _indexToVector;

   // shared globally
   std::vector<THREE_VECTOR> _centers;
   std::vector<double> _radii;
   std::vector<double> _bias;

   // local only
   std::vector<AnatomyCell>& _cells;
   std::vector<std::vector<int> > _nbrDomains;
   
};

#endif
