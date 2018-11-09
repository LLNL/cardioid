#include "sensorFactory.hh"
#include <iostream>
#include <cassert>
#include <string>

#include "object_cc.hh"
#include "Sensor.hh"
#include "PointListSensor.hh"
#include "ActivationTimeSensor.hh"
#include "ActivationAndRecoverySensor.hh"
#include "MinMaxSensor.hh"
#include "DataVoronoiCoarsening.hh"
#include "GradientVoronoiCoarsening.hh"
#include "CaAverageSensor.hh"
#include "MaxDVSensor.hh"
#include "DVThreshSensor.hh"
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
   Sensor* scanActivationAndRecoverySensor(OBJECT* obj, const SensorParms& sp, const Anatomy& anatomy,
                                           const PotentialData&);
   Sensor* scanVoronoiCoarseningSensor(OBJECT* obj, const SensorParms& sp, const Anatomy& anatomy,
                                       const PotentialData&, const Simulate& sim);
   Sensor* scanCaSensor(OBJECT* obj, const SensorParms& sp, const Anatomy& anatomy,
                        const ReactionManager&, const Simulate& sim);
   Sensor* scanMaxDvSensor(OBJECT* obj, const SensorParms& sp, const Anatomy& anatomy,
                           const PotentialData& vdata);
   Sensor* scanDVThreshSensor(OBJECT* obj, SensorParms& sp, const Anatomy& anatomy,
                           const PotentialData& vdata);
   Sensor* scanStateVariableSensor(OBJECT* obj, const SensorParms& sp, const Simulate& sim);
   Sensor* scanECGSensor(OBJECT* obj, const SensorParms& sp, const Simulate& sim);
}


/*!
  @page obj_SENSOR SENSOR object

  Used to write values of interest to a file.

  The following keywords are recognized by all SENSORS.  Additional
  information about each SENSOR, including method specific keywords, can
  be found in the corresponding method subsection.

  @beginkeywords

    @kw{endTime, Time at which the SENSOR stops providing data.
      The eval and print functions will not be called after the end
      time., 1e100 milliseconds}
    @kw{evalRate, Rate in time steps at which the sensor eval function is
      called, 1}
    @kw{method, Choose from
      "activationTime"\,
      "activationAndRecovery"\,
      "averageCa"\,
      "dataVoronoiCoarsening"\,
      "ECG"\,
      "gradientVoronoiCoarsening"\.
      "maxDV"\,
      "DVThresh"\,
      "MinMax"\,
      "pointList"\,
      "stateVariable"\,
      and
      "voronoiCoarsening"
      ,
      No default}
    @kw{printRate}{Rate in time steps at which the sensor print
    function is called.}{1}
    @kw{startTime, Time at which the SENSOR starts providing data.
    The eval and print functions will not be called before the start
    time., -1e100 milliseconds}
    @endkeywords

    @subpage SENSOR_activationTime

    @subpage SENSOR_averageCa

    @subpage SENSOR_dataVoronoiCoarsening
    
    @subpage SENSOR_ECG
    
    @subpage SENSOR_gradientVoronoiCoarsening
    
    @subpage SENSOR_maxDV
    
    @subpage SENSOR_DVThresh
    
    @subpage SENSOR_MinMax
    
    @subpage SENSOR_pointList
    
    @subpage SENSOR_stateVariable
    
    @subpage SENSOR_voronoiCoarsening
*/
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
  else if (method == "activationTime")
     return scanActivationTimeSensor(obj, sp, sim.anatomy_,sim.vdata_);
  else if (method == "activationAndRecovery")
     return scanActivationAndRecoverySensor(obj, sp, sim.anatomy_,sim.vdata_);
  else if (method == "averageCa")
     return scanCaSensor(obj, sp, sim.anatomy_,*sim.reaction_, sim);
#ifdef USE_CUDA
  else if (method == "ECG")
     return scanECGSensor(obj, sp, sim);
#endif
  else if (method == "maxDV")
     return scanMaxDvSensor(obj, sp, sim.anatomy_,sim.vdata_);
  else if (method == "DVThresh")
     return scanDVThreshSensor(obj, sp, sim.anatomy_,sim.vdata_);
  else if (method == "minmax" || method == "MinMax")
     return scanMinMaxSensor(obj, sp, sim.anatomy_,sim.vdata_);
  else if (method == "pointList")
     return scanPointListSensor(obj, sp, sim.anatomy_,sim.vdata_);
  else if (method == "stateVariable")
     return scanStateVariableSensor(obj, sp, sim);
  else if (method == "voronoiCoarsening" 
        || method == "dataVoronoiCoarsening" 
        || method == "gradientVoronoiCoarsening" )
     return scanVoronoiCoarseningSensor(obj, sp, sim.anatomy_,sim.vdata_, sim);


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
   /*!
     @page SENSOR_activationTime SENSOR activationTime method

     Writes a file containing the time at which each cell is first
     activated.  A cell is considered activated when the membrane
     voltage is greater than zero.  Cells that have not been activated
     have an activation time of zero.

     @beginkeywords
     @kw{filename, Name of the pio file into which data is written,
         activationTime}
     @kw{nFiles, The number of physical files for each pio file.
         If nFiles is set to zero cardioid will choose a default value
         that is scaled to the number of MPI ranks in the job.,
         0}
     @endkeywords
    */
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
   /*!
     @page SENSOR_activationAndRecovery SENSOR activationandRecovery method

     Writes a file containing the times at which each cell activates and
     recovers.  Unlike the activation time sensor that only writes the
     first activation, this sensor writes all activation and all
     recovery times since the last time the print function was called.

     A cell is considered activated when the membrane voltage crosses
     the threshhold with positive slope.  Recovery is crossing the
     threshhold with a negative slope.  This sensor is intended to
     satisfy the requirements of the Second N-version Cardiac
     Electrophysiology Benchmark Specification.

     Note that this senor does not write a file with the header usually
     associated with pio files.  The file is formatted to satisfy the
     benchmark spec.  Namely, each record lists the x, y, and z
     coordinates of a cell and all of the activation and recovery times
     for that cell.
     x0 y0 z0 tA1 tR1 tA2 tR2 ...
     x1 y1 z1 tA1 tR1 tA2 tR2 ...
     

     @beginkeywords
     @kw{filename, Name of the pio file into which data is written,
         activationTime}
     @kw{nFiles, The number of physical files for each pio file.
         If nFiles is set to zero cardioid will choose a default value
         that is scaled to the number of MPI ranks in the job.,
         0}
     @kw{threshhold, The voltage above which a cell is considered
     activated., -40 mV}
     @endkeywords
    */
   Sensor* scanActivationAndRecoverySensor(
      OBJECT* obj, const SensorParms& sp, const Anatomy& anatomy,
      const PotentialData& vdata)
   {
      ActivationAndRecoverySensorParms p;
      objectGet(obj, "filename",   p.filename,   "arTime");
      objectGet(obj, "nFiles",     p.nFiles,     "0");
      objectGet(obj, "threshhold", p.threshhold, "-40", "voltage");
      return new ActivationAndRecoverySensor(sp, p, anatomy, vdata);
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

      double maxDistance;
      objectGet(obj, "maxDistance",  maxDistance,  "100000.0");

      string format;
      objectGet(obj, "format", format, "ascii");
      assert( format.compare("ascii")==0 || format.compare("bin")==0 );

      string method;
      objectGet(obj, "method", method, "undefined");
      if ( method == "voronoiCoarsening" ||
           method == "dataVoronoiCoarsening" )
         return new DataVoronoiCoarsening(sp, filename, nFiles, anatomy, cellVec, vdata, sim.commTable_, maxDistance);
      else if( method == "gradientVoronoiCoarsening" )
      {
         string algo;
         objectGet(obj, "algorithm", algo, "comm");
         assert( algo.compare("comm")==0 || algo.compare("nocomm")==0 );
         const bool use_communication_avoiding_algorithm = ( algo.compare("nocomm")==0 );

         return new GradientVoronoiCoarsening(sp, filename, nFiles, anatomy, cellVec, vdata, sim.commTable_, format, maxDistance,
                                              use_communication_avoiding_algorithm);
      }
      assert(false);
      return 0;
   }
}

namespace
{
   Sensor* scanCaSensor(OBJECT* obj, const SensorParms& sp, const Anatomy& anatomy,
                        const ReactionManager& reaction,
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
   /*!
     @page SENSOR_maxDV SENSOR maxDV method

     Prints the minimum and maximum value of $dV_m/dt$ over all active
     cells.  The evalRate is ignored for this sensor.

     @beginkeywords
     @kw(filename, Name for output file., stdout}
     @endkeywords
     
    */
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
   /*!
     @page SENSOR_DVThresh SENSOR DVThresh method

     Stops the simulation if the absolute value of both the minimum and maximum 
     of $dV_m/dt$ over all active cells is less than the given threshold.

     @beginkeywords
     @kw(filename, Name for output file., stdout}
     @endkeywords
     
    */
   Sensor* scanDVThreshSensor(OBJECT* obj, SensorParms& sp, const Anatomy& anatomy,
                           const PotentialData& vdata)
   {
      //get threshold value from object.data
      double val;
      objectGet(obj, "threshold",  val,  "-1.0");
      sp.value = val;
      return new DVThreshSensor(sp, anatomy, vdata, MPI_COMM_WORLD);
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
      objectGet(obj, "allFields", p.allFields, "0");
      objectGet(obj, "cells",    p.cells);
      objectGet(obj, "allCells", p.allCells, "0");
      string outputType; objectGet(obj, "outputType", outputType, "ascii");
      p.binaryOutput =  (outputType != "ascii");

      return new StateVariableSensor(sp, p, sim);      
   }
}
#ifdef USE_CUDA
namespace
{
   Sensor* scanECGSensor(OBJECT* obj, const SensorParms& sp, const Simulate& sim)
   {
      ECGSensorParms p;
      objectGet(obj, "filename",      p.filename,      "ecgData");
      //objectGet(obj, "stencilSize",   p.stencilSize,   "4");
      //objectGet(obj, "nSensorPoints", p.nSensorPoints, "4");
      objectGet(obj, "nFiles",        p.nFiles,        "0");
      objectGet(obj, "kconst",        p.kconst,        "0.8");
      std::vector<std::string> ecgNames;
      objectGet(obj, "ecgPoints",     ecgNames);
      
      p.nSensorPoints=ecgNames.size();
      const int dim=3;
      p.ecgPoints.resize(p.nSensorPoints*dim, 0.0);
      p.ecgNames = ecgNames;
      
      for(int i=0; i<ecgNames.size(); ++i){
          OBJECT* ecgObj = objectFind(ecgNames[i], "POINT");
          objectGet(ecgObj, "x",      p.ecgPoints[i*dim],        "0.0");
          objectGet(ecgObj, "y",      p.ecgPoints[i*dim+1],      "0.0");
          objectGet(ecgObj, "z",      p.ecgPoints[i*dim+2],      "0.0");
          
      }
      return new ECGSensor(sp, p, sim);      
   }
}
#endif
