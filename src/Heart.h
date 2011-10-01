#ifndef HEART_H
#define HEART_H
#include "EPhysModel.h"
#include "MechModel.h"

struct Heart
{
    EPhysModel* ep;   // abstract base class for electrophysiology model
    MechModel* mech;  // abstract base class for mechanical model
};
#endif
