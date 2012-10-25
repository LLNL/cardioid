#include "GradientVoronoiCoarsening.hh"
#include "PerformanceTimers.hh"
#include "pio.h"
#include "ioUtils.h"
#include "Simulate.hh"
#include "CommTable.hh"
using namespace PerformanceTimers;

#include <iostream>
#include <sstream>
#include <iomanip>
#include <list>
#include <string.h>
using namespace std;

// algorithm inspired by reference:
// C.D. Correa, R. Hero, K.-W. Ma, 
// "A comparison of gradient estimation methods for volume
//  rendering on unstructured meshes"
// (see section 3.2)

namespace
{
   /** Concatenates the strings in vv into a single string with a space
    * separating each element.  No space is added to the beginning or
    * end. */
   string concat(const vector<string> vv)
   {
      if (vv.size() == 0)
         return "";
      string tmp = vv[0];
      for (unsigned ii=1; ii<vv.size(); ++ii)
         tmp += " " + vv[ii];
      return tmp;
   }
}

/////////////////////////////////////////////////////////////////////

static const double tol_det=1.e-8;

inline double det3(const double a,const double b,const double c,
            const double d,const double e,const double f,
            const double g,const double h,const double i)
{
   return (a*e*i)+(b*f*g)+(c*d*h)-(c*e*g)-(a*f*h)-(b*d*i);
}

// Cramer's rule
int GradientVoronoiCoarsening::solve3x3(const double s[6], const double r[3], double x[3], const int color)
{
#if 0
   double det1=det3(s[0],s[1],s[2],s[1],s[3],s[4],s[2],s[4],s[5]);
   if( fabs(det1)<tol_det )
   {
      cout<<"WARNING: Bad condition number for 3x3 system: det="<<det1<<", color="<<color<<endl;
      cout<<"A=["<<s[0]<<" "<<s[1]<<" "<<s[2]<<";"
                 <<s[1]<<" "<<s[3]<<" "<<s[4]<<";"
                 <<s[2]<<" "<<s[4]<<" "<<s[5]<<"]"<<endl;
      //cout<<"r=["<<r[0]<<" "<<r[1]<<" "<<r[2]<<"]"<<endl;
      return 1;
   }
   const double det1i=1./det1;
#else
   //assert( invDetMat_.find(color)!=invDetMat_.end() );
   const double det1i=invDetMat_[color];
#endif
   x[0]=det3(r[0],s[1],s[2],r[1],s[3],s[4],r[2],s[4],s[5])*det1i;
   x[1]=det3(s[0],r[0],s[2],s[1],r[1],s[4],s[2],r[2],s[5])*det1i;
   x[2]=det3(s[0],s[1],r[0],s[1],s[3],r[1],s[2],s[4],r[2])*det1i;
   
   return 0;
}

/////////////////////////////////////////////////////////////////////

GradientVoronoiCoarsening::GradientVoronoiCoarsening(const SensorParms& sp,
                                     string filename,
                                     const Anatomy& anatomy,
                                     const vector<Long64>& gid,
                                     const PotentialData& vdata,
                                     const CommTable* commtable,
                                     const string format,
                                     const double max_distance,
                                     const bool use_communication_avoiding_algorithm)
   :Sensor(sp),
    coarsening_(anatomy,gid,commtable),
    filename_(filename),
    anatomy_(anatomy),
    vdata_(vdata),
    comm_(commtable->_comm),
    format_(format),
    max_distance_(max_distance),
    use_communication_avoiding_algorithm_(use_communication_avoiding_algorithm)
{
   const int nLocal = anatomy_.nLocal();
   
   // color local cells
   int ret=coarsening_.bruteForceColoring(max_distance);
   assert( ret>=0 );
   
   dx_.resize(nLocal);
   dy_.resize(nLocal);
   dz_.resize(nLocal);

   coarsening_.colorDisplacements(dx_,dy_,dz_);

   coarsening_.computeRemoteTasks();
   
   const double maxd2=max_distance*max_distance;
   colored_cells_.reserve(27);
   
   
   //int count=0;
   for(int ic=0;ic<nLocal;++ic)
   {
      const int color=coarsening_.getColor(ic);
      if( color>=0 )
      {
         assert( dx_[ic]*dx_[ic]+dy_[ic]*dy_[ic]+dz_[ic]*dz_[ic]<maxd2 ); 
         colored_cells_.push_back(ic);
         //count++;
      }
   }
   
   //cout<<count<<" colored cells"<<endl;
}

void GradientVoronoiCoarsening::computeColorCenterValues(const VectorDouble32& val)
{
   startTimer(sensorCompColorCenterTimer);

   static bool first_time=true;

   // calculate local sums
   
   static list<int> center_indexes;
   static list<int> no_center_indexes;
   static map<int,int> ncellsofcolor;
   
   if( first_time )
   {
      const int nLocal = anatomy_.nLocal();

      for(int ic=0;ic<nLocal;++ic)
      {
         const int color=coarsening_.getColor(ic);
         if( color>=0 )
         {
            const double norm2=(dx_[ic]*dx_[ic]
                               +dy_[ic]*dy_[ic]
                               +dz_[ic]*dz_[ic]); 
            if( norm2<1.e-6 ){
               center_indexes.push_back(ic);
            }
            ncellsofcolor[color]++;
         }
      }
   }
   
   // set sum of n values to 0. for all colors known locally
   if( !use_communication_avoiding_algorithm_ || first_time )
   {
      map<int,int>::const_iterator p=ncellsofcolor.begin();
      while(p!=ncellsofcolor.end())
      {
         const int color=p->first;
      
         valcolors_.setSum(color,ncellsofcolor[color],0.);
         ++p;
      }
   }
   
   // set sum of n values to val[ic] ( (n-1)*0. + val[ic] )
   // for all colors centered locally
   list<int>::const_iterator pi=center_indexes.begin();
   while(pi!=center_indexes.end())
   {
      const int ic=*pi;

      const int color=coarsening_.getColor(ic);
      assert( color>=0 );
      
      valcolors_.setSum(color,ncellsofcolor[color],val[ic]);
      
      ++pi;
   }

   if( !use_communication_avoiding_algorithm_ )
      coarsening_.exchangeAndSum(valcolors_);
   
   first_time=false;

   stopTimer(sensorCompColorCenterTimer);
}

// setup matrix of least square system dX^T W^2 dX grad V = dX^T W^2 dF
void GradientVoronoiCoarsening::setupLSmatrix()
{
   valMat00_.clear();
   valMat01_.clear();
   valMat02_.clear();
   valMat11_.clear();
   valMat12_.clear();
   valMat22_.clear();

   const int nLocalcoloredCells = colored_cells_.size();
   const double minnorm2=1.e-8;

   for(int icc=0;icc<nLocalcoloredCells;++icc)
   {
      int ic=colored_cells_[icc];
      
      const double norm2=(dx_[ic]*dx_[ic]+dy_[ic]*dy_[ic]+dz_[ic]*dz_[ic]); 
      if( norm2>minnorm2) // skip point at distance 0 
      {
         const int color=coarsening_.getColor(ic);
         assert( color>=0 );
         
         const double norm2i=1./norm2; // weighting
         
         valMat00_.add1value(color,dx_[ic]*dx_[ic]*norm2i);
         valMat01_.add1value(color,dx_[ic]*dy_[ic]*norm2i);
         valMat02_.add1value(color,dx_[ic]*dz_[ic]*norm2i);
         valMat11_.add1value(color,dy_[ic]*dy_[ic]*norm2i);
         valMat12_.add1value(color,dy_[ic]*dz_[ic]*norm2i);
         valMat22_.add1value(color,dz_[ic]*dz_[ic]*norm2i);
      }else{
         const int color=coarsening_.getColor(ic);
         assert( color>=0 );
         valMat00_.add1value(color,0.);
         valMat01_.add1value(color,0.);
         valMat02_.add1value(color,0.);
         valMat11_.add1value(color,0.);
         valMat12_.add1value(color,0.);
         valMat22_.add1value(color,0.);
      }
   }
   
   // consolidate matrices over MPI tasks
   vector<LocalSums*> valcolors;
   valcolors.push_back(&valMat00_);
   valcolors.push_back(&valMat01_);
   valcolors.push_back(&valMat02_);
   valcolors.push_back(&valMat11_);
   valcolors.push_back(&valMat12_);
   valcolors.push_back(&valMat22_);
   
   coarsening_.exchangeAndSum(valcolors);
}
   
// setup r.h.s. of least square system dX^T W^2 dX grad V = dX^T W^2 dF
void GradientVoronoiCoarsening::setupLSsystem(const VectorDouble32& val)
{
   startTimer(sensorSetupLSTimer);
   const int nLocalcoloredCells = colored_cells_.size();
   const double minnorm2=1.e-8;
   
   static map<int,double*> bf0;

   static bool first_time = true;
   if( first_time )
   {
      setupLSmatrix();
      
      if( use_communication_avoiding_algorithm_ )
      {
         valRHS0_.clear();
         valRHS1_.clear();
         valRHS2_.clear();
         for(int icc=0;icc<nLocalcoloredCells;++icc)
         {
            int ic=colored_cells_[icc];
      
            const double norm2=(dx_[ic]*dx_[ic]+dy_[ic]*dy_[ic]+dz_[ic]*dz_[ic]); 
            if( norm2>minnorm2) // skip point at distance 0 
            {
               const int color=coarsening_.getColor(ic);
               assert( color>=0 );
         
               const double norm2i=1./norm2; // weighting
         
               valRHS0_.add1value(color,dx_[ic]*norm2i);
               valRHS1_.add1value(color,dy_[ic]*norm2i);
               valRHS2_.add1value(color,dz_[ic]*norm2i);
            }else{
               const int color=coarsening_.getColor(ic);
               assert( color>=0 );
               valRHS0_.add1value(color,0.);
               valRHS1_.add1value(color,0.);
               valRHS2_.add1value(color,0.);
            }
         }
      
         vector<LocalSums*> valcolors;
         valcolors.push_back(&valRHS0_);
         valcolors.push_back(&valRHS1_);
         valcolors.push_back(&valRHS2_);
   
         coarsening_.exchangeAndSum(valcolors);
      
         const std::set<int>& compute_colors(coarsening_.getOwnedColors());
         for(set<int>::const_iterator it = compute_colors.begin();
                                      it!= compute_colors.end();
                                    ++it)
         {
            const int color=(*it);
            double* array=new double[3];
            array[0]=valRHS0_.value(color);
            array[1]=valRHS1_.value(color);
            array[2]=valRHS2_.value(color);
            bf0.insert(pair<int,double*>(color,array));
         }
      }

      first_time=false;
   }

   valRHS0_.clear();
   valRHS1_.clear();
   valRHS2_.clear();

   if( use_communication_avoiding_algorithm_ )
   {
      for(int icc=0;icc<nLocalcoloredCells;++icc)
      {
         int ic=colored_cells_[icc];
      
         const double norm2=(dx_[ic]*dx_[ic]+dy_[ic]*dy_[ic]+dz_[ic]*dz_[ic]); 
         if( norm2>minnorm2) // skip point at distance 0 
         {
            const int color=coarsening_.getColor(ic);
            assert( color>=0 );
         
            const double norm2i=1./norm2; // weighting
         
            const double v=norm2i*val[ic];
            valRHS0_.add1value(color,dx_[ic]*v);
            valRHS1_.add1value(color,dy_[ic]*v);
            valRHS2_.add1value(color,dz_[ic]*v);
         }
      }
      const std::set<int>& compute_colors(coarsening_.getOwnedColors());
      for(set<int>::const_iterator it = compute_colors.begin();
                                   it!= compute_colors.end();
                                 ++it)
      {
         const int color=(*it);
         const double v=-1.*valcolors_.value(color);
         
         valRHS0_.increaseValue(color,v*bf0[color][0]);
         valRHS1_.increaseValue(color,v*bf0[color][1]);
         valRHS2_.increaseValue(color,v*bf0[color][2]);
      }
      
   }else{
      for(int icc=0;icc<nLocalcoloredCells;++icc)
      {
         int ic=colored_cells_[icc];
      
         const double norm2=(dx_[ic]*dx_[ic]+dy_[ic]*dy_[ic]+dz_[ic]*dz_[ic]); 
         if( norm2>minnorm2) // skip point at distance 0 
         {
            const int color=coarsening_.getColor(ic);
            assert( color>=0 );
         
            const double norm2i=1./norm2; // weighting
         
            const double v=norm2i*(val[ic]-valcolors_.value(color));
            valRHS0_.add1value(color,dx_[ic]*v);
            valRHS1_.add1value(color,dy_[ic]*v);
            valRHS2_.add1value(color,dz_[ic]*v);
         }
      }

      // consolidate vector over MPI tasks
      vector<LocalSums*> valcolors;
      valcolors.push_back(&valRHS0_);
      valcolors.push_back(&valRHS1_);
      valcolors.push_back(&valRHS2_);
   
      coarsening_.exchangeAndSum(valcolors);   
   }
   stopTimer(sensorSetupLSTimer);
}

void GradientVoronoiCoarsening::prologComputeLeastSquareGradients()
{
   int myRank;
   MPI_Comm_rank(comm_, &myRank);
   
   if( myRank==0 )
      cout<<"GradientVoronoiCoarsening::prologComputeLeastSquareGradients()..."<<endl;

   nb_excluded_pts_=0;
   
   const std::set<int>& compute_colors = 
      use_communication_avoiding_algorithm_ ? coarsening_.getLocalColors()
                                            : coarsening_.getOwnedColors();

   Long64 nSnapSubLoc = compute_colors.size();
   Long64 nSnapSub;
   MPI_Reduce(&nSnapSubLoc, &nSnapSub, 1, MPI_LONG_LONG, MPI_SUM, 0, comm_);
   
   if( myRank==0 )
      cout<<"nSnapSub="<<nSnapSub<<endl;

   matLS_.clear();
   invDetMat_.clear();
   
   for(set<int>::const_iterator it = compute_colors.begin();
                                it!= compute_colors.end();
                              ++it)
   {
      const int color=(*it);
   
      const int ncells=valMat00_.nValues(color);
      //assert( ncells==valMat01_.nValues(color) );
      //assert( ncells==valMat02_.nValues(color) );
      //assert( ncells==valMat11_.nValues(color) );
      //assert( ncells==valMat12_.nValues(color) );
      //assert( ncells==valMat22_.nValues(color) );
   
      if( ncells>3 ){
         double* a=new double[6];
         a[0]=valMat00_.value(color);
         a[1]=valMat01_.value(color);
         a[2]=valMat02_.value(color);
         a[3]=valMat11_.value(color);
         a[4]=valMat12_.value(color);
         a[5]=valMat22_.value(color);
      
         double det1=det3(a[0],a[1],a[2],a[1],a[3],a[4],a[2],a[4],a[5]);
         if( fabs(det1)>tol_det )
         {
            matLS_.insert(pair<int,double*>(color,a));
         
            invDetMat_.insert(pair<int,double>(color,1./det1));
            
            included_eval_colors_.insert(color);
            
            gradients_[color].reserve(3*printRate()/evalRate());
         }else{
            cout<<"WARNING: unable to compute gradient because of bad condition number of matrix: color "
                <<color<<" will be skipped..."<<endl;
            delete[] a;
         }
         
      }else{
         cout<<"GradientVoronoiCoarsening --- WARNING: exclude "
               "point from coarsened data on task "<<myRank
             <<" (surrounded by "<<ncells<<" cells only)"
             <<endl;
      }
   } 
   
   int npts=((int)compute_colors.size()-(int)included_eval_colors_.size());
   if( npts>0 )cout<<"GradientVoronoiCoarsening --- WARNING: exclude "
                   <<npts<<" points from coarsened data on task "<<myRank
                   <<endl;
   
   MPI_Reduce(&npts, &nb_excluded_pts_, 1, MPI_INT, MPI_SUM, 0, comm_);
   nb_sampling_pts_=nSnapSub-(Long64)nb_excluded_pts_;

   if( myRank==0 )
      cout<<"GradientVoronoiCoarsening: prolog done... nb_sampling_pts_="<<nb_sampling_pts_<<endl;
}

void GradientVoronoiCoarsening::computeLeastSquareGradients(const double current_time,
                                                            const int current_loop)
{
   startTimer(sensorComputeLSTimer);
   static bool first_time = true;
   
   if( first_time )
   {
      prologComputeLeastSquareGradients();
      
      first_time = false;
   }

   for(set<int>::const_iterator it = included_eval_colors_.begin();
                                it!= included_eval_colors_.end();
                              ++it)
   {
      const int color=(*it);
      //cout<<"computeLeastSquareGradients, color="<<color<<endl;
      
      vector<float>& color_gradient(gradients_[color]);
      
      // use first 2 records (6 fields) to store x,y,z and nb. values used for averaging
      if( color_gradient.empty() )
      if( format_.compare("ascii")!=0 )
      {
         const int halfNx = anatomy_.nx()/2;
         const int halfNy = anatomy_.ny()/2;
         const int halfNz = anatomy_.nz()/2;

         const Vector& v = coarsening_.center(color);
         float ix = float(v.x()) - float(halfNx);
         float iy = float(v.y()) - float(halfNy);
         float iz = float(v.z()) - float(halfNz);
      
         color_gradient.push_back(ix);
         color_gradient.push_back(iy);
         color_gradient.push_back(iz);
         
         color_gradient.push_back( float(valcolors_.nValues(color)) );
         // fill with dummies to get a record of length 3
         float dummy=0.;
         color_gradient.push_back(dummy);
         color_gradient.push_back(dummy);
      }
      
      double b[3]={valRHS0_.value(color),
                   valRHS1_.value(color),
                   valRHS2_.value(color)};
      double norm2b=b[0]*b[0]+b[1]*b[1]+b[2]*b[2];          
      
      double g[3]={0.,0.,0.};             
      if( norm2b>1.e-15){
         assert( matLS_.find(color)!=matLS_.end() );
         solve3x3(matLS_[color],b,g,color);
      }
      
      color_gradient.push_back(float(g[0]));
      color_gradient.push_back(float(g[1]));
      color_gradient.push_back(float(g[2]));
   }

#ifdef DEBUG  
   int myRank;
   MPI_Comm_rank(comm_, &myRank);
   if( myRank==0 )
      cout<<"done with GradientVoronoiCoarsening::computeLeastSquareGradients()..."<<endl;
#endif
   stopTimer(sensorComputeLSTimer);
}
   
void GradientVoronoiCoarsening::writeGradients(const string& filename,
                                               const double current_time,
                                               const int current_loop)const
{
   int myRank;
   MPI_Comm_rank(comm_, &myRank);

   PFILE* file = Popen(filename.c_str(), "w", comm_);


   if (myRank == 0)
   {      
      // write header
      int nfiles;
      Pget(file,"ngroup",&nfiles);
      Pprintf(file, "gradV FILEHEADER {\n");
      if( format_.compare("ascii")!=0 )
      {
         Pprintf(file, "  datatype = FIXRECORDBINARY;\n");
         const int lrec    = 3*sizeof(float);
         Pprintf(file, "  lrec = %d;\n", lrec);
         Pprintf(file, "  nfields = %d;\n", 3);
         Pprintf(file, "  nrecords = %llu;\n", nb_sampling_pts_*(times_.size()+2));
         Pprintf(file, "  field_names = gx gy gz\n");
         Pprintf(file, "  field_types = f4 f4 f4;\n");
      }else{
         Pprintf(file, "  datatype = FIXRECORDASCII;\n");
         const int lrec    = 25+13*3*(int)times_.size();
         const int nfields = 4+3*(int)times_.size(); 
         Pprintf(file, "  lrec = %d;\n", lrec);
         Pprintf(file, "  nfields = %d;\n", nfields);
         Pprintf(file, "  nrecords = %llu;\n", nb_sampling_pts_);
         string fieldNames="rx ry rz nvals " + concat(vector<string>(times_.size(), "gx gy gz"));
         Pprintf(file, "  field_names = %s;\n", fieldNames.c_str());
         string fieldTypes="d d d d " + concat(vector<string>(3*times_.size(), "f"));
         Pprintf(file, "  field_types = %s;\n", fieldTypes.c_str());
      }
      Pprintf(file, "  nfiles = %u;\n", nfiles);
      Pprintf(file, "  time = %f; loop = %u;\n", times_[0], current_loop);
      if( times_.size()>1 )
         Pprintf(file, "  nsteps = %d; dt = %f\n", times_.size(), times_[1]-times_[0]);
      Pprintf(file, "  h = %4u  0    0\n", anatomy_.nx());
      Pprintf(file, "        0    %4u  0\n", anatomy_.ny());
      Pprintf(file, "        0    0    %4u;\n", anatomy_.nz());
      Pprintf(file, "}\n\n");
   }
   
   if( format_.compare("ascii")!=0 )
   {
      for(map< int, vector<float> >::const_iterator itn =gradients_.begin();
                                                    itn!=gradients_.end();
                                                  ++itn)
      {
         const int color=(itn->first);

         const short m=coarsening_.multiplicity(color); // usually, m=1
         for(short ii=0;ii<m;ii++)
            Pwrite(&itn->second[0], sizeof(float), itn->second.size(), file);
      }
   }else{
      const int halfNx = anatomy_.nx()/2;
      const int halfNy = anatomy_.ny()/2;
      const int halfNz = anatomy_.nz()/2;
   
      for(set<int>::const_iterator it = included_eval_colors_.begin();
                                   it!= included_eval_colors_.end();
                                 ++it)
      {
         const int color=(*it);

         const Vector& v = coarsening_.center(color);
         int ix = int(v.x()) - halfNx;
         int iy = int(v.y()) - halfNy;
         int iz = int(v.z()) - halfNz;

         const map< int, vector<float> >::const_iterator itn=gradients_.find(color);
         const vector<float>& color_gradient(itn->second);
      
         stringstream ss;
         ss << setw(5)<< right << ix<<" ";
         ss << setw(5)<< right << iy<<" ";
         ss << setw(5)<< right << iz<<" ";
         ss << setw(7)<< right << valcolors_.nValues(color);
      
         ss << setprecision(8);
         for(int it=0;it<times_.size();++it)
         {
            ss <<" "<< setw(12)<< right << color_gradient[3*it+0];
            ss <<" "<< setw(12)<< right << color_gradient[3*it+1];
            ss <<" "<< setw(12)<< right << color_gradient[3*it+2];

         }
         ss << endl;
         string line(ss.str());
         const short m=coarsening_.multiplicity(color);
         for(short ii=0;ii<m;ii++)
            Pwrite(line.c_str(), line.size(), 1, file);
      }
   }
      
   Pclose(file);
}

void GradientVoronoiCoarsening::eval(double time, int loop)
{
   startTimer(sensorEvalTimer);
   
   times_.push_back(time);
   
   computeColorCenterValues(vdata_.VmArray_);
   setupLSsystem(vdata_.VmArray_);
   computeLeastSquareGradients(time, loop);   
   
   stopTimer(sensorEvalTimer);
}

void GradientVoronoiCoarsening::print(double time, int loop)
{
   startTimer(sensorPrintTimer);
   
   int myRank;
   MPI_Comm_rank(comm_, &myRank);

   stringstream name;
   name << "snapshot."<<setfill('0')<<setw(12)<<loop;
   string fullname = name.str();
   if (myRank == 0)
      DirTestCreate(fullname.c_str());
   fullname += "/" + filename_;
   
   writeGradients(fullname,time, loop);   

   times_.clear();
   for(map<int,std::vector<float> >::iterator itg =gradients_.begin();
                                              itg!=gradients_.end();
                                            ++itg)
   {
      (itg->second).clear();
      (itg->second).reserve(3*printRate()/evalRate());
   }
   
   stopTimer(sensorPrintTimer);
}

