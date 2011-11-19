#ifndef SENSOR_FACTORY
#define SENSOR_FACTORY

#include <string>
class Sensor;
class Anatomy;

Sensor* sensorFactory(const std::string& name, const Anatomy& anatomy);

#endif
