#ifndef CONTROL_H
#define CONTROL_H

#include <string>
#include <vector>

struct Control
{
    std::string ephysmodel;
    std::string mechmodel;
    double dt;
    double stimAmplitude;
    int stimLength;
    int stimInterval;
    int savestep;
    int laststep;
    std::vector<int> ngrid_ep;
    std::vector<int> npegrid_ep;
    std::vector<int> ngrid_mech;
    std::vector<int> npegrid_mech;
};
#endif
