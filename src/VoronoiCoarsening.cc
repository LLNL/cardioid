#include <iostream>

#include "VoronoiCoarsening.hh"
#include "pio.h"

VoronoiCoarsening::VoronoiCoarsening(const Anatomy& anatomy,
                                     const std::vector<Long64>& gid,
                                     MPI_Comm comm)
   :anatomy_(anatomy),
    cells_(anatomy.cellArray()),
    comm_(comm),
    indexToVector_(anatomy.nx(), anatomy.ny(), anatomy.nz())
{
   assert( gid.size()>0 );
   
   // initialize centers_ with vectors corresponding to gids passed in
   for(std::vector<Long64>::const_iterator igid =gid.begin();
                                           igid!=gid.end();
                                           ++igid)
   {
      centers_.push_back(indexToVector_(*igid));
   }

   colors_.resize(anatomy.nLocal()); // color only local cells
  
   // color local cells
   int ret=bruteForceColoring();
   assert( ret>=0 );
   
   computeRemoteTasks();
}

// set values of color_ according to index of closest centers
int VoronoiCoarsening::bruteForceColoring()
{
   assert( colors_.size()>0 );
   assert( centers_.size()>0 );
   assert( colors_.size()>0 );
   
   ncolors_.clear();
   local_colors_.clear();
   
   // color one cell at a time
   for (int icell=0; icell<colors_.size(); ++icell)
   {
      double r2Min = 1e30;
      int color = -1;
      Vector r = indexToVector_(cells_[icell].gid_);
      for (int icenter=0; icenter<centers_.size(); ++icenter)
      {
         Vector rij = r - centers_[icenter];
         double r2 = dot(rij, rij);
         if (r2 < r2Min)
         {
            r2Min = r2;
            color = icenter;
         }
      }
      if (color < 0 ){
         std::cerr << "Fail to assign color to cell "<<icell<<std::endl;
         return -1;
      }else{
         colors_[icell]=color;
         ncolors_[color]++;
         local_colors_.insert(color);
      }
   }
   return 0;
}

void VoronoiCoarsening::computeRemoteTasks()
{
   for(std::map<int,int>::iterator itr=ncolors_.begin();itr!=ncolors_.end();++itr)
      assert( local_colors_.find(itr->first)!=local_colors_.end() );

   int nTasks, myRank;
   MPI_Comm_size(comm_, &nTasks);
   MPI_Comm_rank(comm_, &myRank);  

   
   int nlocalcolors=(int)(local_colors_.size());
   int max_nlocalcolors;
   MPI_Allreduce(&nlocalcolors, &max_nlocalcolors, 1, MPI_INT, MPI_MAX, comm_);
   
   // copy local_colors_ into vector
   std::vector<int> local_colors(local_colors_.size());
   copy(local_colors_.begin(),local_colors_.end(),local_colors.begin());
   local_colors.resize(max_nlocalcolors,-1);
   
   // calculate all the remote tasks which share colors with local task
   // (simplest algorithm... may not scale...)
   {
      remote_tasks_.clear();
      std::vector<int> remote_colors(max_nlocalcolors*nTasks,-1);
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
   std::vector<int> nlocal_colors;
   for(std::map<int,int>::iterator itr=ncolors_.begin();itr!=ncolors_.end();++itr)
   {
      nlocal_colors.push_back((int)itr->first);
      nlocal_colors.push_back((int)itr->second);
   }
   nlocal_colors.resize(2*max_nlocalcolors,-1);

   int it=0;
   std::vector< std::vector<int> > nremote_colors_for_task;
   nremote_colors_for_task.resize(remote_tasks_.size());
   for(std::set<int>::const_iterator  itp =remote_tasks_.begin();
                                      itp!=remote_tasks_.end();
                                      ++itp)
   {
      // recv
      nremote_colors_for_task[it].resize(2*max_nlocalcolors);
      MPI_Irecv(&nremote_colors_for_task[it][0], 2*max_nlocalcolors, MPI_INT, *itp, tag, comm_, recvReq+it);

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
   
   // remove -1 at end of nremote_colors_for_task
   for(int ii=0;ii<nremote_colors_for_task.size();ii++)
   {
      std::vector<int>& nremote_colors=nremote_colors_for_task[ii];
      assert( nremote_colors.size()%2==0 );
      while( !nremote_colors.empty() )
      {
         std::vector<int>::iterator ir=nremote_colors.end();
         ir--;
         if( (*ir)==-1 ){
            nremote_colors.pop_back();
            nremote_colors.pop_back();
         }else{
            break;
         }
      }
      it++;
   }
   
   //for(std::set<int>::const_iterator  itp =remote_tasks_.begin();
   //                                   itp!=remote_tasks_.end();
   //                                   ++itp)
   //   std::cout<<"PE "<<myRank<<" exchange data with task "<<*itp<<std::endl;
   
   // just remove local colors that are centered somewhere else
   int tid=0;
   for(std::set<int>::const_iterator  itp =remote_tasks_.begin();
                                      itp!=remote_tasks_.end();
                                      ++itp)
   {
      std::vector<int>& nremote_colors=nremote_colors_for_task[tid];
      int count=0;
      // reduce local_colors_ so that each color is on one task only
      // (the one where the corresponding color is the most present)
      for(int i=0;i<nremote_colors.size()/2;++i)
      {
         const int color=nremote_colors[2*i];
         const int nc   =nremote_colors[2*i+1];
         std::set<int>::const_iterator iti=local_colors_.find(color);
         if( iti!=local_colors_.end() )
         {
            bool should_recv_data=true;
            if( nc>=ncolors_[color] )
            {
               if( nc>ncolors_[color] || myRank>(*itp))
               {
                  local_colors_.erase(iti);
               }
            }
         }
      }
      
      tid++;
   }

   // removed unused communications
   dst_tasks_=remote_tasks_;
   src_tasks_=remote_tasks_;
   
   ncolors_to_send_.clear();
   ncolors_to_recv_.clear();
   
   tid=0;
   for(std::set<int>::const_iterator  itp =remote_tasks_.begin();
                                      itp!=remote_tasks_.end();
                                      ++itp)
   {
      std::vector<int>& nremote_colors=nremote_colors_for_task[tid];
      int count=0;
      // reduce local_colors_ so that each color is on one task only
      // (the one where the corresponding color is the most present)
      for(int i=0;i<nremote_colors.size()/2;++i)
      {
         const int color=nremote_colors[2*i];
         const int nc   =nremote_colors[2*i+1];
         std::set<int>::const_iterator iti=local_colors_.find(color);
         if( iti!=local_colors_.end() )
         {
            bool should_recv_data=true;
            if( nc>=ncolors_[color] )
            {
               if( nc>ncolors_[color] || myRank>(*itp))
               {
                  should_recv_data=false;
               }
            }
            if( should_recv_data )ncolors_to_recv_[*itp]++;
            else                  ncolors_to_send_[*itp]++;
            count++;
         }
      }
      //if(count==0)std::cout<<"PE "<<myRank<<" has no data to exchange with PE "<<*itp<<std::endl;
      
      // remove unnecessay communications
      if( ncolors_to_recv_.find(*itp)==ncolors_to_recv_.end() )src_tasks_.erase(*itp);
      if( ncolors_to_send_.find(*itp)==ncolors_to_send_.end() )dst_tasks_.erase(*itp);
      
      tid++;
   }
#if 0
   for(std::set<int>::const_iterator  itp =src_tasks_.begin();
                                      itp!=src_tasks_.end();
                                      ++itp)
      std::cout<<"PE "<<myRank<<" should receive data from task "<<*itp<<std::endl;
   
   for(std::set<int>::const_iterator  itp =dst_tasks_.begin();
                                      itp!=dst_tasks_.end();
                                      ++itp)
     std::cout<<"PE "<<myRank<<" should send data to task "<<*itp<<std::endl;
#endif
MPI_Barrier(comm_);
 
}

void VoronoiCoarsening::computeColorAverages(const std::vector<double>& val)
{
   // calculate local sums
   valcolors_.clear();
   const int nLocal = colors_.size();
   for(int ic=0;ic<nLocal;++ic)
   {
      valcolors_.add1value(colors_[ic],val[ic]);
   }
   
   // set up send buffer
   PackedData* localpackeddata=new PackedData[ncolors_.size()];
   valcolors_.packData(localpackeddata,ncolors_.size());
   
   // set up buffer to recv remote data
   int ndata2recv=0;
   for(std::map<int,int>::const_iterator p = ncolors_to_recv_.begin();
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
   for(std::set<int>::iterator  itp =src_tasks_.begin();
                                itp!=src_tasks_.end();
                                ++itp)
   {
      // recv
      MPI_Irecv(&remotepackeddata[offset], ncolors_to_recv_[*itp]*sizeof(PackedData), MPI_BYTE, *itp, tag, comm_, recvReq+it);
      offset+=ncolors_to_recv_[*itp];
      it++;
   }
   
   it=0;
   for(std::set<int>::iterator  itp =dst_tasks_.begin();
                                itp!=dst_tasks_.end();
                                ++itp)
   {
      // send
      MPI_Isend(&localpackeddata[0],  ncolors_to_send_[*itp]*sizeof(PackedData), MPI_BYTE, *itp, tag, comm_, sendReq+it);
      it++;
   }

   // wait
   if(recvReq!=0)MPI_Waitall(dst_tasks_.size(), sendReq, MPI_STATUS_IGNORE);
   if(sendReq!=0)MPI_Waitall(src_tasks_.size(), recvReq, MPI_STATUS_IGNORE);
   
   if(recvReq!=0)delete[] recvReq;
   if(sendReq!=0)delete[] sendReq;
   
   // accumulate data in valcolors_
   offset=0;
   for(std::set<int>::iterator  itp =src_tasks_.begin();
                                itp!=src_tasks_.end();
                                ++itp)
   {
      for(int i=0;i<ncolors_to_recv_[*itp];i++)
      {
         int ii=offset+i;
         int color=remotepackeddata[ii].color;
         if( local_colors_.find(color)!=local_colors_.end() )
         {
            valcolors_.addnvalues(color,remotepackeddata[ii].value,remotepackeddata[ii].n);
         }
      }
      
      offset+=ncolors_to_recv_[*itp];
   }
   
   if( ndata2recv>0 )delete[] remotepackeddata;
}


void VoronoiCoarsening::writeAverages(const std::string& filename,
                                      const double current_time,
                                      const int current_loop)const
{
   int myRank;
   MPI_Comm_rank(comm_, &myRank);

   PFILE* file = Popen(filename.c_str(), "w", comm_);

   char fmt[] = "%5d %5d %5d %5d %18.12f";
   int lrec = 55;
   int nfields = 5; 

   Long64 nSnapSub = -1;
   Long64 nSnapSubLoc = local_colors_.size();
   MPI_Allreduce(&nSnapSubLoc, &nSnapSub, 1, MPI_LONG_LONG, MPI_SUM, comm_);

   if (myRank == 0)
   {
      // write header
      int nfiles;
      Pget(file,"ngroup",&nfiles);
      Pprintf(file, "cellViz FILEHEADER {\n");
      Pprintf(file, "  lrec = %d;\n", lrec);
      Pprintf(file, "  datatype = FIXRECORDASCII;\n");
      Pprintf(file, "  nrecords = %llu;\n", nSnapSub);
      Pprintf(file, "  nfields = %d;\n", nfields);
      Pprintf(file, "  field_names = rx ry rz nvals avgVm;\n");
      Pprintf(file, "  field_types = u u u u f;\n" );
      Pprintf(file, "  nfiles = %u;\n", nfiles);
      Pprintf(file, "  time = %f; loop = %u;\n", current_time, current_loop);
      Pprintf(file, "  h = %4u  0    0\n", anatomy_.nx());
      Pprintf(file, "        0    %4u  0\n", anatomy_.ny());
      Pprintf(file, "        0    0    %4u;\n", anatomy_.nz());
      Pprintf(file, "}\n\n");
   }
   
   char line[lrec+1];
   const int halfNx = anatomy_.nx()/2;
   const int halfNy = anatomy_.ny()/2;
   const int halfNz = anatomy_.nz()/2;
   
   for(std::set<int>::const_iterator it = local_colors_.begin();
                                     it!= local_colors_.end();
                                     ++it)
   {
      const int color=(*it);
      const Vector& v = centers_[color];
      int ix = int(v.x()) - halfNx;
      int iy = int(v.y()) - halfNy;
      int iz = int(v.z()) - halfNz;
      
      int l = snprintf(line, lrec, fmt,
                       ix, iy, iz,
                       valcolors_.nValues(color),
                       valcolors_.averageValue(color));
      
      if (myRank == 0 && l>=lrec ){
         std::cerr<<"ERROR: printed record truncated in file "<<filename<<std::endl;
         std::cerr<<"This could be caused by out of range values"<<std::endl;
         std::cerr<<"record="<<line<<std::endl;
         break;
      }
      for (; l < lrec - 1; l++) line[l] = (char)' ';
      line[l++] = (char)'\n';
      assert (l==lrec);
      Pwrite(line, lrec, 1, file);
   }
   
   Pclose(file);
}
