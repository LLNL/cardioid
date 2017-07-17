#ifndef GRID_ASSIGNMENT_OBJECT_H
#define GRID_ASSIGNMENT_OBJECT_H

#include "three_algebra.h"
#include "intQueue.h"

#ifdef __cplusplus
extern "C" 
{
#endif

typedef struct GaoCell_st
{
   int first;
   int burned;
}
GAO_CELL;

/** The GRID_ASSIGNMENT_OBJECT computes the closest center to a given
 * particle coordinate using a grid/flood approach.  The intent of this
 * code is to provide an initial assignment method that scales only as
 * the number of particles to assign.  (I.e., it is independent of the
 * number of centers).
 *
 * To vastly simplify the code we completely ignore periodic boundary
 * conditions.  We can get away with this because the initial assignment
 * doesn't have to be perfect, it only needs to be close.  If we can get
 * a particle into a domain that is close to its correct Voronoi domain
 * then the regular assignment will do the right thing.  */

typedef struct GridAssignmentObject_st
{
   int nCenters;
   const void* centerP;
   int stride;
   int* nextCenter;
   THREE_VECTOR corner;
   int nx, ny, nz;
   double dx, dy, dz;
   GAO_CELL* grid;
   INT_QUEUE* floodQueue;
   INT_QUEUE* burnList;
}
GRID_ASSIGNMENT_OBJECT;


GRID_ASSIGNMENT_OBJECT* gao_init(int nCenters,
				 const void* centerP,
				 int stride);

void gao_destroy(GRID_ASSIGNMENT_OBJECT* gao);

int gao_nearestCenter(GRID_ASSIGNMENT_OBJECT* gao, const THREE_VECTOR r);


#ifdef __cplusplus
}
#endif

#endif
