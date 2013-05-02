#include "sensorFactory.hh"
#include <iostream>
#include <cassert>
#include <string>

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
#include "ECGSensor.hh"
#include "Simulate.hh"
#include "readCellList.hh"

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
                                       const PotentialData&, const Simulate& sim);
   Sensor* scanCaSensor(OBJECT* obj, const SensorParms& sp, const Anatomy& anatomy,
                        const Reaction&, const Simulate& sim);
   Sensor* scanMaxDvSensor(OBJECT* obj, const SensorParms& sp, const Anatomy& anatomy,
                           const PotentialData& vdata);
   Sensor* scanStateVariableSensor(OBJECT* obj, const SensorParms& sp, const Simulate& sim);
   Sensor* scanECGSensor(OBJECT* obj, const SensorParms& sp, const Simulate& sim);
}


Sensor* sensorFactory(const std::string& name, const Simulate& sim)
{
  OBJECT* obj = objectFind(name, "SENSOR");
  string method;
  objectGet(obj, "method", method, "undefined");
  SensorParms sp;
  objectGet(obj, "startTime", sp.startTime, "-1e100", "t"); //time units
  objectGet(obj, "endTime",   sp.endTime,   "1e100", "t");
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
     return scanCaSensor(obj, sp, sim.anatomy_,*sim.reaction_, sim);
  else if (method == "maxDV")
     return scanMaxDvSensor(obj, sp, sim.anatomy_,sim.vdata_);
  else if (method == "voronoiCoarsening" 
        || method == "dataVoronoiCoarsening" 
        || method == "gradientVoronoiCoarsening" )
     return scanVoronoiCoarseningSensor(obj, sp, sim.anatomy_,sim.vdata_, sim);
  else if (method == "stateVariable")
     return scanStateVariableSensor(obj, sp, sim);
  else if (method == "ECG")
     return scanECGSensor(obj, sp, sim);


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
      objectGet(obj, "filename",  p.filename, "activationTime");
      objectGet(obj, "nFiles",    p.nFiles,   "0");
      return new ActivationTimeSensor(sp, p, anatomy, vdata);
   }
}

namespace
{
   Sensor* scanVoronoiCoarseningSensor(OBJECT* obj, const SensorParms& sp, const Anatomy& anatomy,
                                       const PotentialData& vdata,
                                       const Simulate& sim)
   {
      string cellListFilename;
      objectGet(obj, "cellList", cellListFilename, "");

      vector<Long64> cellVec;
      readCellList(cellListFilename, cellVec);

      string filename;
      objectGet(obj, "filename",  filename,  "coarsened_anatomy");

      unsigned nFiles;  objectGet(obj, "nFiles", nFiles, "0");

      double maxdistance;
      objectGet(obj, "maxDistance",  maxdistance,  "100000.0");

      string format;
      objectGet(obj, "format", format, "ascii");
      assert( format.compare("ascii")==0 || format.compare("bin")==0 );

      string method;
      objectGet(obj, "method", method, "undefined");
      if ( method == "voronoiCoarsening" ||
           method == "dataVoronoiCoarsening" )
         return new DataVoronoiCoarsening(sp, filename, nFiles, anatomy, cellVec, vdata, sim.commTable_, maxdistance);
      else if( method == "gradientVoronoiCoarsening" )
      {
         string algo;
         objectGet(obj, "algorithm", algo, "comm");
         assert( algo.compare("comm")==0 || algo.compare("nocomm")==0 );
         const bool use_communication_avoiding_algorithm = ( algo.compare("nocomm")==0 );

         return new GradientVoronoiCoarsening(sp, filename, nFiles, anatomy, cellVec, vdata, sim.commTable_, format, maxdistance,
                                              use_communication_avoiding_algorithm);
      }
      assert(false);
      return 0;
   }
}

namespace
{
   Sensor* scanCaSensor(OBJECT* obj, const SensorParms& sp, const Anatomy& anatomy,
                        const Reaction& reaction,
                        const Simulate& sim)
   {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

      string cellListFilename;
      objectGet(obj, "cellList", cellListFilename, "");

      vector<Long64> cellVec;
      readCellList(cellListFilename, cellVec);

      string filename;
      objectGet(obj, "filename",  filename,  "coarsened_Ca");

      unsigned nFiles; objectGet(obj, "nFiles", nFiles, "0");
      
      return new CaAverageSensor(sp, filename, nFiles, anatomy, cellVec, reaction, sim.commTable_);
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
      assert(object_testforkeyword(obj, "gid") == 0);  // catch obsolete usage 
      objectGet(obj, "nFiles",   p.nFiles,   "0");
      objectGet(obj, "filename", p.filename, "variables");
      objectGet(obj, "cellList", p.cellListFilename, "");
      objectGet(obj, "radius",   p.radius,   "0.0");
      objectGet(obj, "fields",   p.fieldList);
      objectGet(obj, "cells",    p.cells);
      string outputType; objectGet(obj, "outputType", outputType, "ascii");
      p.binaryOutput =  (outputType != "ascii");

      return new StateVariableSensor(sp, p, sim);      
   }
}
namespace
{
   Sensor* scanECGSensor(OBJECT* obj, const SensorParms& sp, const Simulate& sim)
   {
      ECGSensorParms p;
      objectGet(obj, "filename",      p.filename,      "ecgData");
      objectGet(obj, "stencilSize",   p.stencilSize,   "4");
      objectGet(obj, "nSensorPoints", p.nSensorPoints, "4");
      objectGet(obj, "nFiles",        p.nFiles,        "0");
      return new ECGSensor(sp, p, sim);      
   }
}

