#include "sensorFactory.hh"
#include <iostream>
#include <cassert>
#include <string>
#include <sstream>

#include "object_cc.hh"
#include "Sensor.hh"
#include "PointListSensor.hh"
#include "ActivationTimeSensor.hh"
#include "DataVoronoiCoarsening.hh"
#include "GradientVoronoiCoarsening.hh"

using namespace std;

namespace
{
   Sensor* scanPointListSensor(OBJECT* obj, const SensorParms& sp, const Anatomy& anatomy);
   Sensor* scanActivationTimeSensor(OBJECT* obj, const SensorParms& sp, const Anatomy& anatomy);
   Sensor* scanVoronoiCoarseningSensor(OBJECT* obj, const SensorParms& sp, const Anatomy& anatomy);
}


Sensor* sensorFactory(const std::string& name, const Anatomy& anatomy)
{
  OBJECT* obj = objectFind(name, "SENSOR");
  string method;
  objectGet(obj, "method", method, "undefined");
  SensorParms sp;
  objectGet(obj, "printRate", sp.printRate, "1");
  objectGet(obj, "evalRate",  sp.evalRate,  "1");


  if (method == "undefined")
    assert(false);
  else if (method == "pointList")
     return scanPointListSensor(obj, sp, anatomy);
  else if (method == "activationTime")
     return scanActivationTimeSensor(obj, sp, anatomy);
  else if (method == "voronoiCoarsening" 
        || method == "dataVoronoiCoarsening" 
        || method == "gradientVoronoiCoarsening" )
     return scanVoronoiCoarseningSensor(obj, sp, anatomy);

  assert(false); // reachable only due to bad input
  return 0;
}


namespace
{
   Sensor* scanPointListSensor(OBJECT* obj, const SensorParms& sp, const Anatomy& anatomy)
   {
      PointListSensorParms p;
      objectGet(obj, "pointList",   p.pointList);
      objectGet(obj, "startTime",   p.startTime,   "0.0",  "t");
      objectGet(obj, "endTime",     p.endTime,     "-1.0", "t");
      objectGet(obj, "printDerivs", p.printDerivs, "0");
      objectGet(obj, "filename",    p.filename,    "cell");
      objectGet(obj, "dirname",     p.dirname,     "sensorData");
      return new PointListSensor(sp, p, anatomy);
   }
}

namespace
{
   Sensor* scanActivationTimeSensor(OBJECT* obj, const SensorParms& sp, const Anatomy& anatomy)
   {
      ActivationTimeSensorParms p;
      objectGet(obj, "filename",  p.filename,  "activationTime");
      return new ActivationTimeSensor(sp, p, anatomy);
   }
}

namespace
{
   Sensor* scanVoronoiCoarseningSensor(OBJECT* obj, const SensorParms& sp, const Anatomy& anatomy)
   {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

      string cellListFilename;
      objectGet(obj, "cellList", cellListFilename, "");

      // try to open file
      int openfail;
      ifstream input;
      if (myRank == 0)
      {
         input.open(cellListFilename.c_str(),ifstream::in);
         openfail = 0;
         if (!input.is_open())
            openfail = 1;
      }
      MPI_Bcast(&openfail, 1, MPI_INT, 0, MPI_COMM_WORLD);
      if (openfail == 1)
      {
         if (myRank == 0)
            cout << "Could not open cell list file " << cellListFilename << endl;
         return false;
      }

      vector<Long64> cellVec;
      int nSubset;
      if (myRank == 0)
      {
         while (!input.eof()) {
            string query;
            if ( !getline ( input, query ) ) {
                break;
            }
            stringstream ss ( query );
            Long64 igid;
            ss >> igid;
            cellVec.push_back(igid);
         }
         nSubset = cellVec.size();
      }   
      MPI_Bcast(&nSubset, 1, MPI_INT, 0, MPI_COMM_WORLD);
      cellVec.resize(nSubset);
      MPI_Bcast(&cellVec[0], nSubset, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

      string filename;
      objectGet(obj, "filename",  filename,  "coarsened_anatomy");

      string method;
      objectGet(obj, "method", method, "undefined");
      if ( method == "voronoiCoarsening" ||
           method == "dataVoronoiCoarsening" )
         return new DataVoronoiCoarsening(sp, filename, anatomy, cellVec, MPI_COMM_WORLD);
      else if( method == "gradientVoronoiCoarsening" )
         return new GradientVoronoiCoarsening(sp, filename, anatomy, cellVec, MPI_COMM_WORLD);
   }
}
