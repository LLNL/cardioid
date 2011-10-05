#ifndef KORADI_TEST_HH
#define KORADI_TEST_HH

#include <vector>
#include "AnatomyReader.hh"
#include "three_algebra.h"
#include "IndexToVector.hh"

class KoradiTest
{
 public:
   
   KoradiTest(AnatomyReader& anatomy,
	      int nCentersPerTask);
   
   void balanceStep();

 private:

   void pickInitialCenters();
   void initialAssignment();
   void computeRadii();

   
   void allGatherCenters();
   void allReduceRadii();



   int _myRank;
   int _nTasks;

   int _nCentersPerTask;
   int _localOffset;

   IndexToVector _indexToVector;
   std::vector<AnatomyCell>& _cells;
   std::vector<THREE_VECTOR> _centers;
   std::vector<double> _radii;
   
};

#endif
