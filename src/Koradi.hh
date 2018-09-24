#ifndef KORADI_HH
#define KORADI_HH

#include <vector>
#include "AnatomyCell.hh"
#include "Vector.hh"
#include "IndexToVector.hh"
#include "IndexToThreeVector.hh"

class Anatomy;

struct KoradiParms
{
   bool verbose;
   int nCentersPerTask;
   int maxVoronoiSteps;
   int maxSteps;
   int outputRate;
   
   double alphaStep;
   double tolerance;
   double nbrDeltaR;
};

class Koradi
{
 public:
   
   Koradi(Anatomy&, const KoradiParms& parms);
   
   
 private:

   void balanceStep();
   void voronoiBalance();
   void distributeCellsEvenly();
   void pickInitialCenters();
   void assignCells();
   void initialAssignment();
   void sendCellsToDestinations();
   void moveCenters();
   void computeRadii();
   void biasAlpha();
   void printStatistics();
   
   void calculateCellDestinations();
   void exchangeCells();
   void bruteForceDistanceCheck();
   void findNbrDomains();
   void computeLoad(std::vector<double>& load);


   bool verbose_;
   int nCentersPerTask_;
   int maxVoronoiSteps_;
   int maxSteps_;
   int reconditionRate_;
   int outputRate_;
   double alphaStep_;
   double tolerance_;
   double nbrDeltaR_;

   int myRank_;
   int nTasks_;
   int localOffset_;

   IndexToVector indexToVector_;
   IndexToThreeVector indexTo3Vector_;

   // shared globally
   std::vector<Vector> centers_;
   std::vector<double> radii_;
   std::vector<double> alpha_;
   std::vector<double> load_;
   double targetLoad_;
   
   // local only
   std::vector<AnatomyCell>& cells_;
   std::vector<std::vector<int> > nbrDomains_;
   
};

#endif
