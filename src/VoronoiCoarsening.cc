#include <iostream>
#include <sstream>
#include <iomanip>
#include <map>
using namespace std;

#include "VoronoiCoarsening.hh"
#include "pio.h"
#include "ioUtils.h"

VoronoiCoarsening::VoronoiCoarsening(const Anatomy& anatomy,
                                     const vector<Long64>& gid,
                                     MPI_Comm comm)
   :anatomy_(anatomy),
    cells_(anatomy.cellArray()),
    comm_(comm),
    indexToVector_(anatomy.nx(), anatomy.ny(), anatomy.nz())
{
   int myRank;
   MPI_Comm_rank(comm_, &myRank);  
   assert( gid.size()>0 );
   
   // initialize centers_ with vectors corresponding to gids passed in
   std::map<Long64,short> multiplicities;
   for(vector<Long64>::const_iterator igid  = gid.begin();
                                      igid != gid.end();
                                    ++igid)
   {
      std::map<Long64,short>::iterator is=multiplicities.find(*igid);
      
      if( is==multiplicities.end() )
      {
         Long64 color=*igid;
         centers_.push_back(indexToVector_(color));
         multiplicities.insert(make_pair(color, 1));
      }else{
         is->second++;
      }
      //if( myRank==0)cout<<"gid="<<*igid<<endl;
   }
   
   multiplicities_.reserve( multiplicities.size() );
   for(std::map<Long64,short>::const_iterator is =multiplicities.begin();
                                              is!=multiplicities.end();
                                              is++)
   {
      multiplicities_.push_back(is->second);
   }

   if( myRank==0)cout<<"VoronoiCoarsening: number of colors = "<<centers_.size()<<endl;

   colors_.resize(anatomy.nLocal()); // color only local cells
   
   assert( centers_.size()==multiplicities_.size() );
}

// set values of color_ according to index of closest centers
int VoronoiCoarsening::bruteForceColoring()
{
   assert( colors_.size()>0 );
   assert( centers_.size()>0 );
   assert( colors_.size()>0 );
   
   int myRank;
   MPI_Comm_rank(comm_, &myRank);  
   if( myRank==0 )
      cout<<"VoronoiCoarsening: brute force coloring..."<<endl;

   ncolors_.clear();
   local_colors_.clear();
   
   const int ncenters=centers_.size();
   
   // get sub-domain mass center
   Vector domain_center(0.,0.,0.);
   for (int icell=0; icell<colors_.size(); ++icell)
   {
      Vector r = indexToVector_(cells_[icell].gid_);
      domain_center+=r;
   }
   domain_center/=(double)colors_.size();
   
   // get sub-domain radius
   double domain_radius=0.;
   for (int icell=0; icell<colors_.size(); ++icell)
   {
      Vector r = indexToVector_(cells_[icell].gid_);
      Vector rij = r - domain_center;
      double r2 = dot(rij, rij);
      if( r2>domain_radius )domain_radius=r2;
   }
   domain_radius=sqrt(domain_radius);
   //cout<<"Domain center: "<<domain_center<<", Domain radius: "<<domain_radius<<endl;

   // get distance from sub-domain center to closest center
   double rcmin = 1.e30;
   for (int icenter=0; icenter<ncenters; ++icenter)
   {
      Vector rij = domain_center - centers_[icenter];
      double r2 = dot(rij, rij);
      if (r2 < rcmin)
      {
         rcmin = r2;
      }
   }
   rcmin=sqrt(rcmin);
   if(rcmin<domain_radius)rcmin=domain_radius;


   // get centers closest to sub-domain
   const double d2min=1.01*(rcmin+domain_radius)*(rcmin+domain_radius);
   map<int,Vector> close_centers;
   for (int icenter=0; icenter<ncenters; ++icenter)
   {
      Vector rij = domain_center - centers_[icenter];
      double r2 = dot(rij, rij);
      
      if( r2<=d2min)close_centers.insert( pair<int,Vector>(icenter,centers_[icenter]) );
   }
   const int nclosecenters=(int)close_centers.size();
   assert( nclosecenters>0 );
   //std::cout<<"nclosecenters="<<nclosecenters<<std::endl;
   
   // color one cell at a time
   for (int icell=0; icell<colors_.size(); ++icell)
   {
      double r2Min = 1e30;
      int color = -1;
      Vector r = indexToVector_(cells_[icell].gid_);
      
      // loop over closest centers
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
      if (color < 0 ){
         cerr << "Failed to assign color to cell "<<icell<<endl;
         return -1;
      }else{
         colors_[icell]=color;
         ncolors_[color]++;
         local_colors_.insert(color);
      }
   }
#if 0
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
#endif
   if( myRank==0 )
      cout<<"VoronoiCoarsening: brute force coloring done..."<<endl;
   
   return 0;
}

void VoronoiCoarsening::colorDisplacements(std::vector<double>& dx,
                                           std::vector<double>& dy,
                                           std::vector<double>& dz)
{
   for (int icell=0; icell<colors_.size(); ++icell)
   {
      Vector r = indexToVector_(cells_[icell].gid_);
      Vector rij = r - centers_[colors_[icell]];
      dx[icell]=rij.x()*anatomy_.dx();
      dy[icell]=rij.y()*anatomy_.dy();
      dz[icell]=rij.z()*anatomy_.dz();
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

void VoronoiCoarsening::setOwnedColors(const map< int, int* >& nremote_colors_from_task, 
                                       const int size_nremote_colors_from_task)
{   
   assert( centers_.size()>0 );
   
   int myRank;
   MPI_Comm_rank(comm_, &myRank);  

   owned_colors_=local_colors_;
   
   // now remove local colors that are centered somewhere else
   for(set<int>::const_iterator  itp =remote_tasks_.begin();
                                 itp!=remote_tasks_.end();
                               ++itp)
   {
      const map< int, int* >::const_iterator itn=nremote_colors_from_task.find(*itp);
      assert( itn!=nremote_colors_from_task.end() );
      const int* const nremote_colors=itn->second;
      int count=0;
      // reduce owned_colors_ so that each color is on one task only
      // (the one where the corresponding color is the most present)
      for(int i=0;i<size_nremote_colors_from_task;++i)
      {
         const int color=nremote_colors[2*i];
         const int nc   =nremote_colors[2*i+1];
         
         assert( color<(int)centers_.size() );
         
         set<int>::const_iterator iti=owned_colors_.find(color);
         if( iti!=owned_colors_.end() )
         {
            assert( color>=0 );
            if( nc>=ncolors_[color] )
            {
               if( nc>ncolors_[color] || myRank<(*itp))
               {
                  owned_colors_.erase(iti);
                  //cout<<"PE "<<myRank<<" remove color "<<*iti<<endl;
               }
            }
         }
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
}

void VoronoiCoarsening::computeRemoteTasks()
{
   int nTasks, myRank;
   MPI_Comm_size(comm_, &nTasks);
   MPI_Comm_rank(comm_, &myRank);  
   if( myRank==0 )
      cout<<"VoronoiCoarsening: compute remote tasks..."<<endl;

   for(map<int,int>::const_iterator itr =ncolors_.begin();
                                    itr!=ncolors_.end();
                                  ++itr)
      assert( local_colors_.find(itr->first)!=local_colors_.end() );
   
   int nlocalcolors=(int)(local_colors_.size());
   int max_nlocalcolors;
   MPI_Allreduce(&nlocalcolors, &max_nlocalcolors, 1, MPI_INT, MPI_MAX, comm_);
   
   // copy local_colors_ into vector
   vector<int> local_colors(local_colors_.size());
   copy(local_colors_.begin(),local_colors_.end(),local_colors.begin());
   local_colors.resize(max_nlocalcolors,-1);
   
   // calculate all the remote tasks which share colors with local task
   // (simplest algorithm... may not scale...)
   {
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
   }
   
   
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

#if 0
   for(set<int>::const_iterator  itp =remote_tasks_.begin();
                                 itp!=remote_tasks_.end();
                               ++itp)
      cout<<"PE "<<myRank<<" exchange data with task "<<*itp<<endl;
   MPI_Barrier(comm_);
   sleep(1);
#endif

   setOwnedColors(nremote_colors_from_task, max_nlocalcolors);

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
}

void VoronoiCoarsening::exchangeAndSum(LocalSums& valcolors)
{
   // set up send buffer
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
   MPI_Barrier(comm_);
   
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

void VoronoiCoarsening::accumulateValues(const vector<double>& val, LocalSums& valcolors)
{
   valcolors.clear();
   const int nLocal = colors_.size();
   for(int ic=0;ic<nLocal;++ic)
   {
      valcolors.add1value(colors_[ic],val[ic]);
   }
}
