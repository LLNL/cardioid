#include "sensorFactory.hh"
#include <iostream>
#include <cassert>
#include <string>
#include <sstream>

#include "object_cc.hh"
#include "Sensor.hh"
#include "PointListSensor.hh"
#include "ActivationTimeSensor.hh"
#include "MinMaxSensor.hh"
#include "DataVoronoiCoarsening.hh"
#include "GradientVoronoiCoarsening.hh"
#include "CaAverageSensor.hh"
#include "MaxDVSensor.hh"
#include "StateVariableSensor.hh"

#include "Simulate.hh"

using namespace std;

namespace
{
   Sensor* scanPointListSensor(OBJECT* obj, const SensorParms& sp, const Anatomy& anatomy,
                               const PotentialData&);
   Sensor* scanMinMaxSensor(OBJECT* obj, const SensorParms& sp, const Anatomy& anatomy,
                               const PotentialData&);
   Sensor* scanActivationTimeSensor(OBJECT* obj, const SensorParms& sp, const Anatomy& anatomy,
                                    const PotentialData&);
   Sensor* scanVoronoiCoarseningSensor(OBJECT* obj, const SensorParms& sp, const Anatomy& anatomy,
                                       const PotentialData&);
   Sensor* scanCaSensor(OBJECT* obj, const SensorParms& sp, const Anatomy& anatomy,
                        const Reaction&);
   Sensor* scanMaxDvSensor(OBJECT* obj, const SensorParms& sp, const Anatomy& anatomy,
                           const PotentialData& vdata);
   Sensor* scanStateVariableSensor(OBJECT* obj, const SensorParms& sp, const Simulate& sim);
}

// Initialize cellVec with gids listed in file filename
int readCelllist(const string filename, vector<Long64>& cellVec)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   // try to open file
   int openfail;
   ifstream input;
   if (myRank == 0)
   {
      input.open(filename.c_str(),ifstream::in);
      openfail = 0;
      if (!input.is_open())
         openfail = 1;
   }
   MPI_Bcast(&openfail, 1, MPI_INT, 0, MPI_COMM_WORLD);
   if (openfail == 1)
   {
      if (myRank == 0)
         cerr << "Could not open cell list file " << filename << endl;
      return -1;
   }

   int nSubset;
   if (myRank == 0)
   {
      while (!input.eof()) {
         string query;
         if ( !getline ( input, query ) ) {
             break;
         }
         istringstream ss ( query );
         Long64 igid;
         while( ss >> igid ){
            assert( igid>=0 );
            cellVec.push_back(igid);
         }
      }
      nSubset = cellVec.size();
   }   
   MPI_Bcast(&nSubset, 1, MPI_INT, 0, MPI_COMM_WORLD);
   cellVec.resize(nSubset);
   MPI_Bcast(&cellVec[0], nSubset, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
   
   return 0;
}

Sensor* sensorFactory(const std::string& name, const Simulate& sim)
{
  OBJECT* obj = objectFind(name, "SENSOR");
  string method;
  objectGet(obj, "method", method, "undefined");
  SensorParms sp;
  objectGet(obj, "printRate", sp.printRate, "1");
  objectGet(obj, "evalRate",  sp.evalRate,  "-1");
  if(sp.evalRate == -1)sp.evalRate=sp.printRate;

  if (method == "undefined")
    assert(false);
  else if (method == "pointList")
     return scanPointListSensor(obj, sp, sim.anatomy_,sim.vdata_);
  else if (method == "minmax" || method == "MinMax")
     return scanMinMaxSensor(obj, sp, sim.anatomy_,sim.vdata_);
  else if (method == "activationTime")
     return scanActivationTimeSensor(obj, sp, sim.anatomy_,sim.vdata_);
  else if (method == "averageCa")
     return scanCaSensor(obj, sp, sim.anatomy_,*sim.reaction_);
  else if (method == "maxDV")
     return scanMaxDvSensor(obj, sp, sim.anatomy_,sim.vdata_);
  else if (method == "voronoiCoarsening" 
        || method == "dataVoronoiCoarsening" 
        || method == "gradientVoronoiCoarsening" )
     return scanVoronoiCoarseningSensor(obj, sp, sim.anatomy_,sim.vdata_);
  else if (method == "stateVariables")
     return scanStateVariableSensor(obj, sp, sim);


   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   if(myRank==0)cerr<<"Sensor ERROR: unknown method "<<method<<endl;

   assert(false); // reachable only due to bad input
   return 0;
}


namespace
{
   Sensor* scanPointListSensor(OBJECT* obj, const SensorParms& sp, const Anatomy& anatomy,
                               const PotentialData& vdata)
   {
      PointListSensorParms p;
      objectGet(obj, "pointList",   p.pointList);
      objectGet(obj, "startTime",   p.startTime,   "0.0",  "t");
      objectGet(obj, "endTime",     p.endTime,     "-1.0", "t");
      objectGet(obj, "printDerivs", p.printDerivs, "0");
      objectGet(obj, "filename",    p.filename,    "cell");
      objectGet(obj, "dirname",     p.dirname,     "sensorData");
      return new PointListSensor(sp, p, anatomy, vdata);
   }
}

namespace
{
   Sensor* scanMinMaxSensor(OBJECT* obj, const SensorParms& sp, const Anatomy& anatomy,
                               const PotentialData& vdata)
   {
      MinMaxSensorParms p;
      objectGet(obj, "startTime",   p.startTime,   "0.0",  "t");
      objectGet(obj, "endTime",     p.endTime,     "-1.0", "t");
      objectGet(obj, "filename",    p.filename,    "minmax");
      objectGet(obj, "dirname",     p.dirname,     "sensorData");
      return new MinMaxSensor(sp, p, anatomy, vdata);
   }
}

namespace
{
   Sensor* scanActivationTimeSensor(OBJECT* obj, const SensorParms& sp, const Anatomy& anatomy,
                                    const PotentialData& vdata)
   {
      ActivationTimeSensorParms p;
      objectGet(obj, "filename",  p.filename,  "activationTime");
      return new ActivationTimeSensor(sp, p, anatomy, vdata);
   }
}

namespace
{
   Sensor* scanVoronoiCoarseningSensor(OBJECT* obj, const SensorParms& sp, const Anatomy& anatomy,
                                       const PotentialData& vdata)
   {
      string cellListFilename;
      objectGet(obj, "cellList", cellListFilename, "");

      vector<Long64> cellVec;
      readCelllist(cellListFilename, cellVec);

      string filename;
      objectGet(obj, "filename",  filename,  "coarsened_anatomy");

      double maxdistance;
      objectGet(obj, "maxDistance",  maxdistance,  "100000.0");

      string format;
      objectGet(obj, "format", format, "ascii");
      assert( format.compare("ascii")==0 || format.compare("bin")==0 );

      string method;
      objectGet(obj, "method", method, "undefined");
      if ( method == "voronoiCoarsening" ||
           method == "dataVoronoiCoarsening" )
         return new DataVoronoiCoarsening(sp, filename, anatomy, cellVec, vdata, MPI_COMM_WORLD, maxdistance);
      else if( method == "gradientVoronoiCoarsening" )
      {
         string algo;
         objectGet(obj, "algorithm", algo, "comm");
         assert( algo.compare("comm")==0 || algo.compare("nocomm")==0 );
         const bool use_communication_avoiding_algorithm = ( algo.compare("nocomm")==0 );

         return new GradientVoronoiCoarsening(sp, filename, anatomy, cellVec, vdata, MPI_COMM_WORLD, format, maxdistance,
                                              use_communication_avoiding_algorithm);
      }
      assert(false);
      return 0;
   }
}

namespace
{
   Sensor* scanCaSensor(OBJECT* obj, const SensorParms& sp, const Anatomy& anatomy,
                        const Reaction& reaction)
   {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

      string cellListFilename;
      objectGet(obj, "cellList", cellListFilename, "");

      vector<Long64> cellVec;
      readCelllist(cellListFilename, cellVec);

      string filename;
      objectGet(obj, "filename",  filename,  "coarsened_Ca");

      return new CaAverageSensor(sp, filename, anatomy, cellVec, reaction, MPI_COMM_WORLD);
   }
}

namespace
{
   Sensor* scanMaxDvSensor(OBJECT* obj, const SensorParms& sp, const Anatomy& anatomy,
                           const PotentialData& vdata)
   {
      string filename;
      objectGet(obj, "filename",  filename,  "cout");

      return new MaxDVSensor(sp, anatomy, vdata, MPI_COMM_WORLD, filename);
   }
}
namespace
{
   Sensor* scanStateVariableSensor(OBJECT* obj, const SensorParms& sp, const Simulate& sim)
   {
      StateVariableSensorParms p;
      objectGet(obj, "gid",   p.gidCenter, "-1");
      assert(p.gidCenter >= 0);
      objectGet(obj, "radius",   p.radius, "0.0");
      objectGet(obj, "fields",   p.fieldList);
      objectGet(obj, "startTime",   p.startTime,   "0.0",  "t");
      objectGet(obj, "endTime",     p.endTime,     "-1.0", "t");
      objectGet(obj, "dirname",     p.dirname,     "stateSensorData");
      return new StateVariableSensor(sp, p, sim);      
   }
}
