#ifndef SENSOR_FACTORY
#define SENSOR_FACTORY

#include <string>
class Sensor;
class Simulate;

Sensor* sensorFactory(const std::string& name, const Simulate& sim);

#endif
