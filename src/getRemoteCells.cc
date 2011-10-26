#include "getRemoteCells.hh"
#include <vector>

#include "Simulate.hh"
#include "Anatomy.hh"
#include "GridRouter.hh"
#include "HaloExchange.hh"

using namespace std;

void getRemoteCells(Simulate& sim, MPI_Comm comm)
{
   Anatomy& anatomy = sim.anatomy_;

   int nx = sim.nx_;
   int ny = sim.ny_;
   int nz = sim.nz_;
   
   vector<Long64> myCells(anatomy.size());
   for (unsigned ii=0; ii<anatomy.size(); ++ii)
      myCells[ii] = anatomy.gid(ii);
   
   sim.router_= new GridRouter(myCells, nx, ny, nz, MPI_COMM_WORLD);
   HaloExchange<AnatomyCell> cellExchange(sim.router_->sendMap(), sim.router_->commTable());

   anatomy.nLocal() = anatomy.size();
   cellExchange.execute(anatomy.cellArray(), anatomy.nLocal());
   anatomy.nRemote() = anatomy.size() - anatomy.nLocal();
}

