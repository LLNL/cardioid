#include "BoxStimulus.hh"
#include "Anatomy.hh"
#include <cmath>
#include <mpi.h>
#include <iostream>

using namespace std;

BoxStimulus::BoxStimulus(const BoxStimulusParms& p, const Anatomy& anatomy, Pulse* pulse, const std::string& name)
: Stimulus(p.baseParms),
  pulse_(pulse)
{
   // Loop through local points and store the indices of any that are in
   // the stimulated region.
   PinnedVector<int> stimList;
   for (unsigned ii=0; ii<anatomy.nLocal(); ++ii)
   {
      Tuple gt = anatomy.globalTuple(ii);
      if (gt.x() > p.xMin && gt.x() < p.xMax &&
          gt.y() > p.yMin && gt.y() < p.yMax &&
          gt.z() > p.zMin && gt.z() < p.zMax )
         stimList.push_back(ii);
   }
   
   // Print number of gids (compute cells) within each box stimulus, this helps you verify that you box stimuli lie completely within tissue
   // For example, if your box is 20x20x20, then you would expect to see 19^3 = 6859 gids inside box stimulus (remember your fenceposts, it is not 20^3 gids, it is 19^3)
   int myrank;								
   int nLocal = stimList.size();						// Local # of box stimulus gids within tissue; the _ at end means this variable only apply to current class, not others (just a naming convention)
   int nGlobal;									// Global # of box stimulus gids within tissue
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);					// Get the current rank (process) #
   MPI_Reduce(&nLocal, &nGlobal, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);	// Add local rank sum to global sum accross all ranks:  [1:  # local gids w/in tissue] [2:  # global gids w/in tissue] [3:  # elements to send from this process] [4:  MPI datatype of each element] [5:  MPI operation to perform] [6:  rank of receiving process w/in communicator] [7:  The MPI communicator]
   if (myrank == 0)								// If at root rank (rank 0), print out global # of cells w/in tissue for this box stimulus
   {
   	cout << "# of tissue compute cells within box stimulus \"" << name << "\" = " << nGlobal << " (# cells should be (lx-1)x(ly-1)x(lz-1) b/c of fencepost condition, where lx, ly, and lz are length of box in discrete points in x, y, and z, respectively)" << endl;
   }
   stimListTransport_.setup(std::move(stimList));
}

int BoxStimulus::subClassStim(double time,
                              Managed<ArrayView<double>> dVmDiffusion)
{
   double value = pulse_->eval(time);
   if (value == 0)
      return 0;
   void addStimulus(OnDevice<ArrayView<double>> dVmDiffusion, OnDevice<ConstArrayView<int>> stimListRaw, const double value);
   addStimulus(dVmDiffusion, stimListTransport_, value);
   return 1;
}

int BoxStimulus::nStim()
{
   ConstArrayView<int> stimList = stimListTransport_.raw();
   return stimList.size();
}
