#include "VoronoiCoarsening.hh"

#include <sstream>
#include <iomanip>
#include <algorithm>

#include "pio.h"
#include "ioUtils.h"
#include "CommTable.hh"
#include "mpiTpl.hh"
#include "GridAssignmentObject.h"
#include "mpiUtils.h"
#include "IndexToThreeVector.hh"

using namespace std;

//#define DEBUG

static void preScreenSensorPoints(vector<Long64>& sensorPoints, const Anatomy& anatomy, MPI_Comm comm);





Vector VoronoiCoarsening::getDomaincenter()const
{
   const int ncells=(int)anatomy_.nLocal();

   // get sub-domain mass center
   Vector domain_center(0.,0.,0.);
   for (int icell=0; icell<ncells; ++icell)
   {
      Vector r = indexToVector_(anatomy_.gid(icell));
      domain_center+=r;
   }
   domain_center/=(double)ncells;
   
   return domain_center;
}      

double VoronoiCoarsening::getDomainRadius(const Vector& domain_center)const
{
   const int ncells=(int)anatomy_.nLocal();

   double domain_radius=0.;
   for (int icell=0; icell<ncells; ++icell)
   {
      Vector r = indexToVector_(anatomy_.gid(icell));
      Vector rij = r - domain_center;
      double r2 = dot(rij, rij);
      if( r2>domain_radius )domain_radius=r2;
   }
   domain_radius=sqrt(domain_radius);
   //cout<<"Domain center: "<<domain_center<<", Domain radius: "<<domain_radius<<endl;

   return domain_radius;
}

map<int,Vector> VoronoiCoarsening::getCloseCenters(const Vector& domain_center, 
                                                   const double d2min)const
{
      map<int,Vector> close_centers;
      for(map<int, Long64>::const_iterator icenter =colorToGidMap_.begin();
                                           icenter!=colorToGidMap_.end();
                                         ++icenter)
      {
         
         Vector ri = indexToVector_(icenter->second);
         Vector rij = domain_center - ri;
         double r2 = dot(rij, rij);
         
         if (r2<=d2min)
            close_centers.insert( make_pair(icenter->first, ri) );
      }
      const int nclosecenters=(int)close_centers.size();
      assert( nclosecenters>0 );
      //std::cout<<"nclosecenters="<<nclosecenters<<std::endl;

   return close_centers;
}


// 1. Screens out any duplicate or non-tissue point in the vector of
//    sensorPoints.
// 2. Assigns a "color" to each sensor point.  Records the assignment in
//    colorToGidMap_.
// 3. Populates ownedColors_ as the set of all colors that correspond
//    to centers that are local on this task.
// 4. Calls gaoForceColoring to
//    4.1 Assign a color to all cells on the task.  (cells that are
//        farther than maxDistance from their center are assigned
//        color = -1.
//    4.2 Build the set localColors_ that contains all colors that
//        are present on the local task (including those colors where
//        the center is remote).
//    4.3 Populate the nColors_ map that gives the number of local cells
//        for each (local) color.
// 5. Sets up comm.
VoronoiCoarsening::VoronoiCoarsening(const Anatomy& anatomy,
                                     vector<Long64>& sensorPoint,
                                     const double maxDistance,
                                     const CommTable* commtable)
: anatomy_(anatomy),
  comm_(commtable->_comm),
  commTable_(commtable),
  indexToVector_(anatomy.nx(), anatomy.ny(), anatomy.nz())
{
   int myRank;
   MPI_Comm_rank(comm_, &myRank);  
   assert( sensorPoint.size() > 0 );

   preScreenSensorPoints(sensorPoint, anatomy, comm_);

   if (myRank==0)
      cout << "VoronoiCoarsening: number of sensor points = "
           << sensorPoint.size() << endl;
   
   // Because the sensor points are now pre-screened, the color to
   // sensor point mapping is now trivial.  I.e., the color of
   // sensorPoint[iColor] is iColor.
   for (unsigned color=0; color<sensorPoint.size(); ++color)
      colorToGidMap_.insert(make_pair(color, sensorPoint[color]));

   
   set<Long64> localCells;
   for (unsigned ii=0; ii<anatomy.nLocal(); ++ii)
      localCells.insert(anatomy.gid(ii));
   for (map<int, Long64>::iterator iter=colorToGidMap_.begin();
        iter!=colorToGidMap_.end(); ++iter)
      if (localCells.count(iter->second) == 1)
         ownedColors_.insert(iter->first);
   
   MPI_Barrier(comm_);
   
   gaoColoring(maxDistance, sensorPoint);
   computeRemoteTasks();   
}

// Initializes:
// - cell_colors_
// - ncolors_
// - localColors_
// (only for cells within maxDistance from a center, other cells take color -1)
int VoronoiCoarsening::gaoColoring(const double maxDistance,
                                   const vector<Long64>& sensorPoint)
{
   timestampBarrier("Starting VoronoiCoarsening::gaoColoring", comm_);

   assert( colorToGidMap_.size()>0 );
   // Check that the mapping from colors to sensorPoint is in fact
   // trivial since we will exploit that fact in gaoColoring.  I think
   // I'd like to get to a day where we don't need colorToGidMap_ and
   // just take advantage of the implicit mapping.  That is a job for
   // another day.
   for (map<int, Long64>::const_iterator iter=colorToGidMap_.begin();
        iter!=colorToGidMap_.end(); ++iter)
      assert(iter->second == sensorPoint[iter->first]);

   int myRank;
   MPI_Comm_rank(comm_, &myRank);  
   IndexToThreeVector indexTo3Vector(anatomy_.nx(), anatomy_.ny(), anatomy_.nz());

   ncolors_.clear();
   localColors_.clear();
   cell_colors_.resize(anatomy_.nLocal(), -1); // color only local cells

   const double r2Max=3.*maxDistance*maxDistance/
      (anatomy_.dx()*anatomy_.dx()
       +anatomy_.dy()*anatomy_.dy()
       +anatomy_.dz()*anatomy_.dz());

   vector<THREE_VECTOR> sensorLocation(sensorPoint.size());
   for (unsigned ii=0; ii<sensorPoint.size(); ++ii)
      sensorLocation[ii] = indexTo3Vector(sensorPoint[ii]);

   GRID_ASSIGNMENT_OBJECT* gao = gao_init(sensorLocation.size(),
                                          (const void*) &sensorLocation[0],
                                          sizeof(THREE_VECTOR));
   
   const int nLocal = anatomy_.nLocal();
   for (unsigned iCell=0; iCell<nLocal; ++iCell)
   {
      THREE_VECTOR rCell = indexTo3Vector(anatomy_.gid(iCell));
      int nearestSensor = gao_nearestCenter(gao, rCell);
      THREE_VECTOR rSensor = sensorLocation[nearestSensor];
      
      double rr = DIFFSQ(rCell, rSensor);
      if (rr > r2Max)
         cell_colors_[iCell] = -1;
      else
      {
         cell_colors_[iCell] = nearestSensor;
         ncolors_[nearestSensor]++;
         localColors_.insert(nearestSensor);
      }
   }

   gao_destroy(gao);
   timestampBarrier("Finished VoronoiCoarsening::gaoColoring", comm_);

   return 0;
}





// Compute the displacement between each colored cell and the
// center cell for that color.
void VoronoiCoarsening::colorDisplacements(std::vector<double>& dx,
                                           std::vector<double>& dy,
                                           std::vector<double>& dz)
{
#ifdef DEBUG
   int myRank;
   MPI_Comm_rank(comm_, &myRank);  
#endif
   const int ncells=(int)anatomy_.nLocal();

   for (int icell=0; icell<ncells; ++icell)
   {
      int color=cell_colors_[icell];
      if(color>=0)
      {
         Vector r = indexToVector_(anatomy_.gid(icell));
         Vector rij = r - indexToVector_(colorToGidMap_[color]);
         dx[icell]=rij.x()*anatomy_.dx();
         dy[icell]=rij.y()*anatomy_.dy();
         dz[icell]=rij.z()*anatomy_.dz();
#ifdef DEBUG
         const double norm2=(dx[icell]*dx[icell]+dy[icell]*dy[icell]+dz[icell]*dz[icell]); 
         if( norm2<1.e-8 )
         {
            assert( ownedColors_.size()>0 );
            
            if( ownedColors_.find(color)==ownedColors_.end() )
            {
               cout<<"colorDisplacements --- ERROR, myRank="<<myRank
                   <<", color="<<color
                   <<", gid="<<anatomy_.gid(icell)
                   <<", r="<<r<<", center="<<colorToGidMap_[color]<<endl;
               for(set<int>::const_iterator it =ownedColors_.begin();
                                            it!=ownedColors_.end();
                                          ++it)
               {
                  cout<<"rank="<<myRank<<" owns color "<<*it<<endl;
               } 
            }
            assert( ownedColors_.find(color)!=ownedColors_.end() );
         }
#endif
      }
   }
}

void VoronoiCoarsening::setupComm(const map< int, int* >& nremote_colors_from_task, 
                                  const int size_nremote_colors_from_task)
{
   int myRank;
   MPI_Comm_rank(comm_, &myRank);  

   dst_tasks_=remote_tasks_;
   src_tasks_=remote_tasks_;
   
   // removed unused communications
   ncolors_to_send_.clear();
   ncolors_to_recv_.clear();
   
   for(set<int>::const_iterator  itp =remote_tasks_.begin();
                                 itp!=remote_tasks_.end();
                               ++itp)
   {
      const map< int, int* >::const_iterator itn=nremote_colors_from_task.find(*itp);
      assert( itn!=nremote_colors_from_task.end() );
      const int* const nremote_colors=itn->second;
      for(int i=0;i<size_nremote_colors_from_task;++i)
      {
         const int color=nremote_colors[2*i];
         const int nc   =nremote_colors[2*i+1];
         set<int>::const_iterator iti=localColors_.find(color);
         
         // do something only if I know that color locally
         //if( iti!=localColors_.end() )
         if( color>=0 )
         {
            ncolors_to_recv_[*itp]++;
         }else{
            break;
         }
      }
      ncolors_to_send_[*itp]=(int)localColors_.size();
   }
}

void VoronoiCoarsening::computeRemoteTasks()
{
   int myRank;
   MPI_Comm_rank(comm_, &myRank);  
   timestampBarrier("Starting VoronoiCoarsening:computeRemoteTasks", comm_);

#if 0
   for(map<int,int>::const_iterator itr =ncolors_.begin();
                                    itr!=ncolors_.end();
                                  ++itr)
      assert( localColors_.find(itr->first)!=localColors_.end() );
#endif
   
   int nlocalcolors=(int)(localColors_.size());
   int max_nlocalcolors;
   MPI_Allreduce(&nlocalcolors, &max_nlocalcolors, 1, MPI_INT, MPI_MAX, comm_);
   if( myRank==0 )
      cout<<"VoronoiCoarsening: max_nlocalcolors/task="<<max_nlocalcolors<<endl;
   
   // copy localColors_ into vector
   vector<int> local_colors(localColors_.size());
   copy(localColors_.begin(),localColors_.end(),local_colors.begin());
   local_colors.resize(max_nlocalcolors,-1);

#if 1
   int nTasks;
   MPI_Comm_size(comm_, &nTasks);
   // calculate all the remote tasks which share colors with local task
   // (simplest algorithm... may not scale...)
   
   remote_tasks_.clear();
   vector<int> remote_colors(max_nlocalcolors*nTasks,-1);
   MPI_Allgather(&local_colors[0],  max_nlocalcolors, MPI_INT, 
                 &remote_colors[0], max_nlocalcolors, MPI_INT, comm_);
   
   // get set of tasks I need to talk to
   for(int pe=0;pe<nTasks;pe++)
   if( pe!=myRank ){
      int* remote_colors_pe=&remote_colors[max_nlocalcolors*pe];
      for(int j=0;j<max_nlocalcolors;j++){
         if( ncolors_.find(remote_colors_pe[j])!=ncolors_.end() ){
            remote_tasks_.insert(pe);
            break;
         }
      }
   }
#else
   
   int nRecvTasks=commTable_->_recvTask.size();
   int nSendTasks=commTable_->_sendTask.size();
   int nTasks=nRecvTasks;
   vector<int> remote_colors(max_nlocalcolors*nRecvTasks,-1);
   
   {
      
      std::vector<MPI_Request> recvReq(nRecvTasks);
      std::vector<MPI_Request> sendReq(nSendTasks);
   
      MPI_Request* precvReq = &recvReq[0];      
      const int tag = 181818;
      for (unsigned ii=0; ii< commTable_->_recvTask.size(); ++ii)
      {
         unsigned sender = commTable_->_recvTask[ii];
         MPI_Irecv(&remote_colors[ii*max_nlocalcolors], max_nlocalcolors, MPI_INT, sender, tag, commTable_->_comm, precvReq+ii);
      }
 
      MPI_Request* psendReq = &sendReq[0];
      for (unsigned ii=0; ii<commTable_->_sendTask.size(); ++ii)
      {
         unsigned target = commTable_->_sendTask[ii];
         MPI_Isend(&local_colors[0], max_nlocalcolors, MPI_INT, target, tag, commTable_->_comm, psendReq+ii);
      }
   }
   // get set of tasks I need to talk to
   for(int ii=0;ii<commTable_->_recvTask.size();ii++){
      unsigned sender = commTable_->_recvTask[ii];
      int* remote_colors_pe=&remote_colors[max_nlocalcolors*ii];
      for(int j=0;j<max_nlocalcolors;j++){
         if( ncolors_.find(remote_colors_pe[j])!=ncolors_.end() ){
            remote_tasks_.insert(sender);
            break;
         }
      }
   }
#endif   
   
   int tag=822;
   MPI_Request* recvReq=0;
   MPI_Request* sendReq=0;
   if( remote_tasks_.size()>0 ){
      recvReq=new MPI_Request[remote_tasks_.size()];
      sendReq=new MPI_Request[remote_tasks_.size()];
   }

   // pack colors and their count into vector
   vector<int> nlocal_colors;
   for(map<int,int>::const_iterator itr =ncolors_.begin();
                                    itr!=ncolors_.end();
                                  ++itr)
   {
      nlocal_colors.push_back((int)itr->first);
      nlocal_colors.push_back((int)itr->second);
   }
   assert( nlocal_colors.size()<=2*max_nlocalcolors );

   nlocal_colors.resize(2*max_nlocalcolors,-1);

   map< int, int* > nremote_colors_from_task;
   int it=0;
   for(set<int>::const_iterator  itp =remote_tasks_.begin();
                                 itp!=remote_tasks_.end();
                               ++itp)
   {
      // recv
      int* buffer=new int[2*max_nlocalcolors];
      MPI_Irecv(&buffer[0], 2*max_nlocalcolors, MPI_INT, *itp, tag, comm_, recvReq+it);
      nremote_colors_from_task.insert(pair<int,int*> (*itp,buffer));

      // send
      MPI_Isend(&nlocal_colors[0],  2*max_nlocalcolors, MPI_INT, *itp, tag, comm_, sendReq+it);
      
      it++;
   }

   // wait
   if( remote_tasks_.size()>0 ){
      MPI_Waitall(remote_tasks_.size(), sendReq, MPI_STATUS_IGNORE);
      MPI_Waitall(remote_tasks_.size(), recvReq, MPI_STATUS_IGNORE);
   
      delete[] recvReq;
      delete[] sendReq;
   }

#ifdef DEBUG
   //for(set<int>::const_iterator  itp =remote_tasks_.begin();
   //                              itp!=remote_tasks_.end();
   //                            ++itp)
   //   cout<<"PE "<<myRank<<" exchange data with task "<<*itp<<endl;
#endif

   setupComm(nremote_colors_from_task, max_nlocalcolors);

   for(set<int>::const_iterator  itp =remote_tasks_.begin();
                                 itp!=remote_tasks_.end();
                               ++itp)
   {
      map< int, int* >::const_iterator itn=nremote_colors_from_task.find(*itp);
      assert( itn!=nremote_colors_from_task.end() );
      delete[] itn->second;
   }
   timestampBarrier("Finished VoronoiCoarsening:computeRemoteTasks", comm_);

   //cout<<"Task "<<myRank<<" receives data from "<<src_tasks_.size()
   //                     <<" tasks and sends data to "<<dst_tasks_.size()<<endl;
}

void VoronoiCoarsening::exchangeAndSum(LocalSums& valcolors)
{
   // set up send buffer
   if( valcolors.size()==0 )
   {
      assert( ncolors_to_recv_.size()==0 );
      assert( src_tasks_.size()==0 );
      assert( dst_tasks_.size()==0 );
      return;
   }

   PackedData* localpackeddata=new PackedData[valcolors.size()];
   valcolors.packData(localpackeddata,valcolors.size());
   
   // set up buffer to recv remote data
   int ndata2recv=0;
   for(map<int,int>::const_iterator p = ncolors_to_recv_.begin();
                                    p!= ncolors_to_recv_.end();
                                  ++p)
      ndata2recv+=p->second;
   
   PackedData* remotepackeddata=0;
   if( ndata2recv>0 )
      remotepackeddata=new PackedData[ndata2recv];
   
   // exchange local sums
   MPI_Request* recvReq=0;
   if(src_tasks_.size()>0)recvReq=new MPI_Request[src_tasks_.size()];
   MPI_Request* sendReq=0;
   if(dst_tasks_.size()>0)sendReq=new MPI_Request[dst_tasks_.size()];
   int tag=823;
   int it=0;
   int offset=0;
   for(set<int>::const_iterator itp =src_tasks_.begin();
                                itp!=src_tasks_.end();
                              ++itp)
   {
      // recv
      MPI_Irecv(&remotepackeddata[offset], 
                ncolors_to_recv_[*itp]*sizeof(PackedData), 
                MPI_BYTE, *itp, tag, comm_, recvReq+it);
      offset+=ncolors_to_recv_[*itp];
      it++;
   }
   
   it=0;
   for(set<int>::const_iterator itp =dst_tasks_.begin();
                                itp!=dst_tasks_.end();
                              ++itp)
   {
      // send
      MPI_Isend(&localpackeddata[0], ncolors_to_send_[*itp]*sizeof(PackedData), MPI_BYTE, *itp, tag, comm_, sendReq+it);
      it++;
   }

   // wait
   if(sendReq!=0)MPI_Waitall(dst_tasks_.size(), sendReq, MPI_STATUS_IGNORE);
   if(recvReq!=0)MPI_Waitall(src_tasks_.size(), recvReq, MPI_STATUS_IGNORE);
   //MPI_Barrier(comm_);
   
   if(recvReq!=0)delete[] recvReq;
   if(sendReq!=0)delete[] sendReq;
   
   // accumulate data in valcolors
   offset=0;
   for(set<int>::const_iterator itp =src_tasks_.begin();
                                itp!=src_tasks_.end();
                              ++itp)
   {
      for(int i=0;i<ncolors_to_recv_[*itp];i++)
      {
         int ii=offset+i;
         int color=remotepackeddata[ii].color;
         if( localColors_.find(color)!=localColors_.end() )
         {
            valcolors.addnvalues(color,remotepackeddata[ii].value,
                                       remotepackeddata[ii].n);
         }
      }
      
      offset+=ncolors_to_recv_[*itp];
   }
   
   delete[] localpackeddata;
   if( ndata2recv>0 )delete[] remotepackeddata;
}

void VoronoiCoarsening::exchangeAndSum(vector<LocalSums*> valcolors)
{
   const unsigned short nvect=(unsigned short)valcolors.size();
   if( valcolors[0]->size()==0 )
   {
      //cout<<"WARNING: valcolors.size()=0 in VoronoiCoarsening::exchangeAndSum()"<<endl;
      assert( ncolors_to_recv_.size()==0 );
      assert( src_tasks_.size()==0 );
      assert( dst_tasks_.size()==0 );
      return;
   }
   
   // set up send buffer
   PackedData* localpackeddata=new PackedData[nvect*valcolors[0]->size()];
   for(unsigned short i=0;i<nvect;i++)
      valcolors[i]->packData(localpackeddata+i*valcolors[i]->size(),valcolors[i]->size());
   
   // set up buffer to recv remote data
   int ndata2recv=0;
   for(map<int,int>::const_iterator p = ncolors_to_recv_.begin();
                                    p!= ncolors_to_recv_.end();
                                  ++p)
      ndata2recv+=p->second;
   
   ndata2recv*=nvect;
   
   PackedData* remotepackeddata=0;
   if( ndata2recv>0 )
      remotepackeddata=new PackedData[ndata2recv];
   
   // exchange local sums
   MPI_Request* recvReq=0;
   if(src_tasks_.size()>0)recvReq=new MPI_Request[src_tasks_.size()];
   MPI_Request* sendReq=0;
   if(dst_tasks_.size()>0)sendReq=new MPI_Request[dst_tasks_.size()];
   int tag=824;
   int it=0;
   int offset=0;
   for(set<int>::const_iterator itp =src_tasks_.begin();
                                itp!=src_tasks_.end();
                              ++itp)
   {
      const int ndata=nvect*ncolors_to_recv_[*itp];
      
      // recv
      MPI_Irecv(&remotepackeddata[offset], 
                ndata*sizeof(PackedData), 
                MPI_BYTE, *itp, tag, comm_, recvReq+it);
      offset+=ndata;
      it++;
   }
   
   it=0;
   for(set<int>::const_iterator itp =dst_tasks_.begin();
                                itp!=dst_tasks_.end();
                              ++itp)
   {
      // send
      MPI_Isend(&localpackeddata[0], nvect*ncolors_to_send_[*itp]*sizeof(PackedData), MPI_BYTE, *itp, tag, comm_, sendReq+it);
      it++;
   }

   // wait
   if(sendReq!=0)MPI_Waitall(dst_tasks_.size(), sendReq, MPI_STATUS_IGNORE);
   if(recvReq!=0)MPI_Waitall(src_tasks_.size(), recvReq, MPI_STATUS_IGNORE);
   //MPI_Barrier(comm_);
   
   if(recvReq!=0)delete[] recvReq;
   if(sendReq!=0)delete[] sendReq;
   
   // accumulate data in valcolors
   offset=0;
   for(set<int>::const_iterator itp =src_tasks_.begin();
                                itp!=src_tasks_.end();
                              ++itp)
   {
      for(unsigned short j=0;j<nvect;j++){
         for(int i=0;i<ncolors_to_recv_[*itp];i++)
         {
            int ii=offset+i;
            int color=remotepackeddata[ii].color;
            if( localColors_.find(color)!=localColors_.end() )
            {
               valcolors[j]->addnvalues(color,remotepackeddata[ii].value,
                                        remotepackeddata[ii].n);
            }
         }
      
         offset+=ncolors_to_recv_[*itp];
      }
   }
   
   delete[] localpackeddata;
   if( ndata2recv>0 )delete[] remotepackeddata;
}

void VoronoiCoarsening::accumulateValues(ro_array_ptr<double> val, LocalSums& valcolors)
{
   valcolors.clear();
   const int nLocal = cell_colors_.size();
   for(int ic=0;ic<nLocal;++ic)
   {
      if( cell_colors_[ic]>=0 )
         valcolors.add1value(cell_colors_[ic],val[ic]);
   }
}

/// Perform checks on the user supplied list of sensorPoints.
/// - Remove any duplicates
/// - Remove and points that do not correspond to tissue cells.
/// On return, sensorPoints will be sorted and free of any duplicates or
/// non-tissue cells.
void preScreenSensorPoints(vector<Long64>& sensorPoints, const Anatomy& anatomy, MPI_Comm comm)
{
   int myRank;
   MPI_Comm_rank(comm, &myRank);  

   // Remove duplicate points.
   sort(sensorPoints.begin(), sensorPoints.end());
   vector<Long64>::iterator
      uniqEnd = unique(sensorPoints.begin(), sensorPoints.end());
   set<Long64> duplicatePoints(uniqEnd, sensorPoints.end());
   sensorPoints.erase(uniqEnd, sensorPoints.end());

   // count appearances of sensorPoints in Anatomy
   set<Long64> localCells;
   for (unsigned ii=0; ii<anatomy.nLocal(); ++ii)
      localCells.insert(anatomy.gid(ii));
   vector<int> count(sensorPoints.size());
   for (unsigned ii=0; ii<sensorPoints.size(); ++ii)
      count[ii] = localCells.count(sensorPoints[ii]);
   allReduce(count, MPI_INT, MPI_SUM, comm);
   
   // erase sensorPoints not in the Anatomy
   vector<Long64> nonTissuePoints;
   for (vector<int>::iterator iter=count.begin(); iter!=count.end();)
   {
      assert(*iter <= 1);
      if (*iter == 0)
      {
         int dist = distance(count.begin(), iter);
         iter = count.erase(iter);
         nonTissuePoints.push_back(*(sensorPoints.begin()+dist));
         sensorPoints.erase(sensorPoints.begin()+dist);
      }
      else
         ++iter;
   }
   
   if (myRank != 0)
      return;
         
   cout << "Duplicate points removed:  " << duplicatePoints.size() << endl;
   cout << "Non-tissue points removed: " << nonTissuePoints.size() << endl;
   
}

