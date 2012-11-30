#include <iostream>
#include <sstream>
#include <iomanip>
#include <map>
using namespace std;

#include "VoronoiCoarsening.hh"
#include "pio.h"
#include "ioUtils.h"
#include "CommTable.hh"

//#define DEBUG

VoronoiCoarsening::VoronoiCoarsening(const Anatomy& anatomy,
                                     const vector<Long64>& gid,
                                     const CommTable* commtable)
   :anatomy_(anatomy),
    comm_(commtable->_comm),
    commTable_(commtable),
    indexToVector_(anatomy.nx(), anatomy.ny(), anatomy.nz())
{
   int myRank;
   MPI_Comm_rank(comm_, &myRank);  
   assert( gid.size()>0 );
   if( myRank==0)cout<<"VoronoiCoarsening: number of gids = "<<gid.size()<<endl;
   
   const unsigned global_ncolors = (int)gid.size();
   vector<short> gidfound(global_ncolors,0);
   const int ncells=(int)anatomy_.nLocal();
   
   owned_colors_.clear();
   for (unsigned color=0; color<global_ncolors; ++color)
   {
      const Long64& rgid = gid[color];
      for (unsigned icell=0; icell<ncells; ++icell)
      {
         if (anatomy_.gid(icell) == rgid)
         {
            gidfound[color] = 1;
            owned_colors_.insert(color);
         }
      }
   }

   vector<short> gidsum(global_ncolors,0);
   MPI_Allreduce(&gidfound[0], &gidsum[0], global_ncolors, MPI_SHORT, MPI_SUM, comm_);
   
   std::map<Long64,int> included_ids;
   for (unsigned color=0; color<global_ncolors; ++color)
   {
      if (gidsum[color] == 0)
      {
         if (myRank == 0)
            cout << "WARNING: VoronoiCoarsening could not find non-zero cell type with gid "
                 << gid[color] << "!  Skipping sensor point. color="<<color << endl;
      }else{
         std::map<Long64,int>::iterator is=included_ids.find(gid[color]);
         if( is==included_ids.end() )
         {
            included_ids.insert(pair<Long64,int>(gid[color],color));
            multiplicities_.insert(make_pair(color, 1));
            
            centers_.insert(pair<int,Vector>(color,indexToVector_(gid[color])));
         }else{
            int cc=is->second;
            if (myRank == 0)
               cout << "WARNING: VoronoiCoarsening, gid "
                    << gid[color] << " already found before... Remove color "<<cc<< endl;
            assert( multiplicities_.find(is->second)!=multiplicities_.end() );
            multiplicities_[cc]++;
            owned_colors_.erase(color);
         }
         //if (myRank == 0)
         //   cout << "gid "<< gid[color] << " is of color "<<color << endl;
      }
   }
#if 0
   for(set<int>::const_iterator  itp =owned_colors_.begin();
                                 itp!=owned_colors_.end();
                               ++itp)
      cout<<"PE "<<myRank<<" owns color "<<*itp<<endl;
   MPI_Barrier(comm_);
   sleep(1);
#endif   
   cout<<"ECG sensors: PE "<<myRank<<" owns "
       <<owned_colors_.size()<<" sensor points "<<endl;

   if( myRank==0)cout<<"VoronoiCoarsening: number of colors = "<<centers_.size()<<endl;

   cell_colors_.resize(anatomy.nLocal(),-1); // color only local cells
   
   if( centers_.size()!=multiplicities_.size() )
   {
      cout<<"centers_.size()="<<centers_.size()
          <<", multiplicities_.size()="<<multiplicities_.size()<<endl;
   }
   
   assert( centers_.size()==multiplicities_.size() );
}

// set values of color_ according to index of closest centers
// (only for cells with max_distance from a center, other cells take color -1)
int VoronoiCoarsening::bruteForceColoring(const double max_distance)
{
   assert( centers_.size()>0 );
   
   int myRank;
   MPI_Comm_rank(comm_, &myRank);  
   if( myRank==0 )
      cout<<"VoronoiCoarsening: brute force coloring..."<<endl;

   ncolors_.clear();
   local_colors_.clear();
   
   const int ncells=(int)anatomy_.nLocal();
   if( ncells>0 )
   {
      if( myRank==0 )cout<<"VoronoiCoarsening: ncenters="<<centers_.size()<<endl;
      
      // get sub-domain mass center
      Vector domain_center(0.,0.,0.);
      for (int icell=0; icell<ncells; ++icell)
      {
         Vector r = indexToVector_(anatomy_.gid(icell));
         domain_center+=r;
      }
      domain_center/=(double)ncells;
      
      // get sub-domain radius
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

      const double r2max=3.*max_distance*max_distance/
                        (anatomy_.dx()*anatomy_.dx()
                        +anatomy_.dy()*anatomy_.dy()
                        +anatomy_.dz()*anatomy_.dz());
      // get distance from sub-domain center to closest center
      double extra_radius = sqrt(r2max);
      if( extra_radius>domain_radius ){
         const double domain_radius2=domain_radius*domain_radius;
         double extra_radius2=extra_radius*extra_radius;
         for(map<int,Vector>::const_iterator icenter =centers_.begin();
                                             icenter!=centers_.end();
                                           ++icenter)
         {
            Vector rij = domain_center - icenter->second;
            double r2 = dot(rij, rij);
            if (r2 < extra_radius2)
            {
               extra_radius2 = r2;
            }
            if( extra_radius2<domain_radius2 )break;
         }
         extra_radius=sqrt(extra_radius2);
         if(extra_radius<domain_radius)extra_radius=domain_radius;
         //std::cout<<"extra_radius="<<extra_radius<<std::endl;
      }

      // get centers closest to sub-domain to reduce cost of coloring later
      const double d2min=1.01*(extra_radius+domain_radius)*(extra_radius+domain_radius);
      map<int,Vector> close_centers;
      for(map<int,Vector>::const_iterator icenter =centers_.begin();
                                          icenter!=centers_.end();
                                        ++icenter)
      {
         Vector rij = domain_center - icenter->second;
         double r2 = dot(rij, rij);
         
         if( r2<=d2min)close_centers.insert( *icenter );
      }
      const int nclosecenters=(int)close_centers.size();
      assert( nclosecenters>0 );
      //std::cout<<"nclosecenters="<<nclosecenters<<std::endl;
      
      //std::cout<<"r2max="<<r2max<<std::endl;
      
      // color one cell at a time
      // color = index of closest sensor point (center)
      for (int icell=0; icell<ncells; ++icell)
      {
         double r2Min = 1e30;
         int color = -1;
         Vector r = indexToVector_(anatomy_.gid(icell));
         
         // loop over closest centers only
         for(map<int,Vector>::const_iterator itr =close_centers.begin();
                                             itr!=close_centers.end();
                                           ++itr)
         {
            Vector rij = r - itr->second;
            double r2 = dot(rij, rij);
            if (r2 < r2Min)
            {
               r2Min = r2;
               color = itr->first;
            }
         }
         
         if( r2Min>r2max ){
            cell_colors_[icell]=-1;
         }else if (color < 0 ){
            cerr << "Failed to assign color to cell "<<icell<<endl;
            return -1;
         }else{
            cell_colors_[icell]=color;
            ncolors_[color]++;
            local_colors_.insert(color);
         }
      }
   }

#ifdef DEBUG
   if( max_distance>99999. )
   {
      int nlocal=0;
      for(map<int,int>::const_iterator itr =ncolors_.begin();
                                       itr!=ncolors_.end();
                                     ++itr)
         nlocal+=itr->second;
      
      int ntotal;
      MPI_Allreduce(&nlocal, &ntotal, 1, MPI_INT, MPI_SUM, comm_);
      if( ntotal!=anatomy_.nGlobal() )
      {
         cerr<<"ERROR in VoronoiCoarsening::bruteForceColoring(): ntotal colors="<<ntotal<<", anatomy_.nGlobal()="<<anatomy_.nGlobal()<<endl;
      }
      assert( ntotal==anatomy_.nGlobal() );
   }
#endif
   if( myRank==0 )
      cout<<"VoronoiCoarsening: brute force coloring done..."<<endl;

   return 0;
}

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
         Vector rij = r - centers_[color];
         dx[icell]=rij.x()*anatomy_.dx();
         dy[icell]=rij.y()*anatomy_.dy();
         dz[icell]=rij.z()*anatomy_.dz();
#ifdef DEBUG
         const double norm2=(dx[icell]*dx[icell]+dy[icell]*dy[icell]+dz[icell]*dz[icell]); 
         if( norm2<1.e-8 )
         {
            if( owned_colors_.size()==0 )
            {
               cout<<"colorDisplacements --- ERROR owned_colors_.size(), myRank="<<myRank
                   <<", color="<<color
                   <<", gid="<<anatomy_.gid(icell)<<endl;
               cout<<"myRank="<<myRank<<", r="<<r<<", center="<<centers_[color]<<endl;
            }
            assert( owned_colors_.size()>0 );
            
            if( owned_colors_.find(color)==owned_colors_.end() )
            {
               cout<<"colorDisplacements --- ERROR, myRank="<<myRank
                   <<", color="<<color
                   <<", gid="<<anatomy_.gid(icell)
                   <<", r="<<r<<", center="<<centers_[color]<<endl;
               for(set<int>::const_iterator it =owned_colors_.begin();
                                            it!=owned_colors_.end();
                                          ++it)
               {
                  cout<<"rank="<<myRank<<" owns color "<<*it<<endl;
               } 
            }
            assert( owned_colors_.find(color)!=owned_colors_.end() );
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
         set<int>::const_iterator iti=local_colors_.find(color);
         
         // do something only if I know that color locally
         //if( iti!=local_colors_.end() )
         if( color>=0 )
         {
            ncolors_to_recv_[*itp]++;
         }else{
            break;
         }
      }
      ncolors_to_send_[*itp]=(int)local_colors_.size();
   }
}

void VoronoiCoarsening::computeRemoteTasks()
{
   int myRank;
   MPI_Comm_rank(comm_, &myRank);  
   if( myRank==0 )
      cout<<"VoronoiCoarsening: compute remote tasks..."<<endl;

#if 0
   for(map<int,int>::const_iterator itr =ncolors_.begin();
                                    itr!=ncolors_.end();
                                  ++itr)
      assert( local_colors_.find(itr->first)!=local_colors_.end() );
#endif
   
   int nlocalcolors=(int)(local_colors_.size());
   int max_nlocalcolors;
   MPI_Allreduce(&nlocalcolors, &max_nlocalcolors, 1, MPI_INT, MPI_MAX, comm_);
   if( myRank==0 )
      cout<<"VoronoiCoarsening: max_nlocalcolors/task="<<max_nlocalcolors<<endl;
   
   // copy local_colors_ into vector
   vector<int> local_colors(local_colors_.size());
   copy(local_colors_.begin(),local_colors_.end(),local_colors.begin());
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
   if( myRank==0 )
      cout<<"VoronoiCoarsening: compute remote tasks done..."<<endl;

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
         if( local_colors_.find(color)!=local_colors_.end() )
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
            if( local_colors_.find(color)!=local_colors_.end() )
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

void VoronoiCoarsening::accumulateValues(const VectorDouble32& val, LocalSums& valcolors)
{
   valcolors.clear();
   const int nLocal = cell_colors_.size();
   for(int ic=0;ic<nLocal;++ic)
   {
      if( cell_colors_[ic]>=0 )
         valcolors.add1value(cell_colors_[ic],val[ic]);
   }
}
