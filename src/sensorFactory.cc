#include "sensorFactory.hh"
#include <iostream>
#include <cassert>
#include "object_cc.hh"
#include "Sensor.hh"
#include "PointListSensor.hh"

using namespace std;

namespace
{
  Sensor* scanPointListSensor(OBJECT* obj, const Anatomy& anatomy);
}


Sensor* sensorFactory(const std::string& name, const Anatomy& anatomy)
{
  OBJECT* obj = objectFind(name, "SENSOR");
  string method;
  objectGet(obj, "method", method, "undefined");
  if (method == "undefined")
    assert(false);
  else if (method == "pointlist")
    return scanPointListSensor(obj, anatomy);
  
  assert(false); // reachable only due to bad input
}


namespace
{
   Sensor* scanPointListSensor(OBJECT* obj, const Anatomy& anatomy)
   {
      PointListSensorParms p;
      objectGet(obj, "pointlist", p.pointlist);
      objectGet(obj, "startTime", p.startTime, "0.0");
      objectGet(obj, "endTime", p.endTime, "-1.0");
      objectGet(obj, "printRate", p.printRate, "1");
      objectGet(obj, "filebase", p.filebase, "sensor.pointlist");
      return new PointListSensor(p, anatomy);
   }
}
