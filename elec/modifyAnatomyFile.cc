#include <string>
#include <sstream>
#include <iomanip>
#include <cassert>

#include "pio.h"
#include "object.h"
#include "BucketOfBits.hh"
#include "object_cc.hh"
#include "readPioFile.hh"
#include "heap.h"

using namespace std;

MPI_Comm COMM_LOCAL = MPI_COMM_WORLD;

string printToWidth(const double value, const int width) {
   stringstream ss;
   int currentPrecision = width;
   for (int currentPrecision=width; currentPrecision>0; --currentPrecision) {
      ss << setw(width) << setprecision(currentPrecision) << value;
      if (ss.str().size() <= width)
      {
         break;
      }
      else
      {
         ss.str("");
         ss.clear();
      }
   }
   return ss.str();
}

//#ifdef BGQ
//#define sym(name) name
//#else
#define sym(name) name ## _
//#endif

extern "C" {
   void sym(dsytrd)(const char* uplo, const int* n, double* A, const int* lda, double* d, double* e, double* tau, double* work, const int* lwork, int* info);
   void sym(dorgtr)(const char* uplo, const int* n, double* A, const int* lda, const double* tau, double* work, const int* lwork, int* info);
   void sym(dsteqr)(const char* compz, const int* n, double* d, double* e, double* Z, const int* ldz, double* work, int* info);
}

void symmetricTensorDecomp(const double inputT[3][3], double vT[3][3], double diag[3])
{
   for (int ii=0; ii<3; ii++)
   {
      for (int jj=0; jj<3; jj++)
      {
         vT[ii][jj] = inputT[ii][jj];
      }
   }

   const char compz='V', uplo='L';
   const int n=3;
   int info;
   double e[n-1], tau[n-1];

   int lwork;
   {
      lwork = -1;
      double work[1];
      sym(dsytrd)(&uplo,&n, &vT[0][0], &n, diag, e, tau, work, &lwork, &info);
      assert(info == 0);
      lwork = work[0];
   }
   {
      double work[lwork];
      sym(dsytrd)(&uplo,&n, &vT[0][0], &n, diag, e, tau, work, &lwork, &info);
      assert(info == 0);
   }
   {
      lwork = -1;
      double work[1];
      sym(dorgtr)(&uplo,&n, &vT[0][0], &n, tau, work, &lwork, &info);
      assert(info == 0);
      lwork = work[0];
   }
   {
      double work[lwork];
      sym(dorgtr)(&uplo,&n, &vT[0][0], &n, tau, work, &lwork, &info);
      assert(info == 0);
   }
   {
      double work[2*n-1];
      sym(dsteqr)(&compz, &n, diag, e, &vT[0][0], &n, work, &info);
      assert(info == 0);
   }
}

class IRegion
{
 public:
   virtual bool contains(const double xx, const double yy, const double zz) const = 0;
   virtual ~IRegion() {}
};

class Box : public IRegion
{
 public:
   virtual bool contains(const double xx, const double yy, const double zz) const
   {
      return (xlow_ <= xx && xx <= xlow_+xextent_)
         &&  (ylow_ <= yy && yy <= ylow_+yextent_)
         &&  (zlow_ <= zz && zz <= zlow_+zextent_)
         ;
   }
   Box(const double xlow,const double ylow,const double zlow,const double xextent, const double yextent, const double zextent) :
   xlow_(xlow),
   ylow_(ylow),
   zlow_(zlow),
   xextent_(xextent),
   yextent_(yextent),
   zextent_(zextent)
   {}
 private:
   double xlow_;
   double ylow_;
   double zlow_;
   double xextent_;
   double yextent_;
   double zextent_;
};

class Ball : public IRegion
{
 public:
   virtual bool contains(const double xx, const double yy, const double zz) const
   {
      return (  (xx-xcenter_)*(xx-xcenter_)
              + (yy-ycenter_)*(yy-ycenter_)
              + (zz-zcenter_)*(zz-zcenter_)
      ) <= radius_*radius_;
   }
   Ball(const double xcenter,const double ycenter, const double zcenter, const double radius) :
   xcenter_(xcenter),
   ycenter_(ycenter),
   zcenter_(zcenter),
   radius_(radius)
   {}
 private:
   double xcenter_;
   double ycenter_;
   double zcenter_;
   double radius_;
};

class EllipsoidInBox : public IRegion
{
 public:
   virtual bool contains(const double xx, const double yy, const double zz) const
   {
      double xrescaled = (xx - (xlow_ + xextent_/2))/xextent_;
      double yrescaled = (yy - (ylow_ + yextent_/2))/yextent_;
      double zrescaled = (zz - (zlow_ + zextent_/2))/zextent_;

      return (  xrescaled*xrescaled
              + yrescaled*yrescaled
              + zrescaled*zrescaled
      ) <= 1;
   }
   EllipsoidInBox(const double xlow,const double ylow,const double zlow,const double xextent, const double yextent, const double zextent) :
   xlow_(xlow),
   ylow_(ylow),
   zlow_(zlow),
   xextent_(xextent),
   yextent_(yextent),
   zextent_(zextent)
   {}
 private:
   double xlow_;
   double ylow_;
   double zlow_;
   double xextent_;
   double yextent_;
   double zextent_;
};

class Edit {
 public:
   IRegion* region;
   int cellType;
   double gil;
   double git;
   double gin;
   bool shouldReplaceCellTypes;
   vector<int> replaceCellTypes;
   Edit() {
      cellType = -1;
      gil = -1;
      git = -1;
      gin = -1;
      shouldReplaceCellTypes = false;
   }
};
      
int main(int argc, char** argv)
{
   int nTasks, myRank;
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   
   if (myRank == 0)
   {
      object_set("files", "modifyAnatomy.data");
      object_compile();
      printf("\nContents of object database:\n"
             "----------------------------------------------------------------------\n");
      object_print_all(stdout);
      printf("----------------------------------------------------------------------\n"
             "End of object database\n\n");
   }
   object_Bcast(0, COMM_LOCAL);

   OBJECT* modifyObj = objectFind("modifyAnatomy", "MODSPEC");
   
   
   //query the object database for things we need to be doing
   string indir;
   objectGet(modifyObj, "indir", indir, "snapshot.initial");
   string outdir;
   objectGet(modifyObj, "outdir", outdir, "outputAnatomy");
   vector<string> editNames;
   objectGet(modifyObj, "edits", editNames);
   int heapSize;
   objectGet(modifyObj, "heap", heapSize, "500");
   heap_start(heapSize);
   
   vector<Edit> edits;
   for (int ii=0; ii<editNames.size(); ii++)
   {
      OBJECT* editObj = objectFind(editNames[ii], "EDIT");
      Edit edit;
      objectGet(editObj, "cellType", edit.cellType, "-1");
      objectGet(editObj, "gil", edit.gil, "-1");
      objectGet(editObj, "git", edit.git, "-1");
      objectGet(editObj, "gin", edit.gin, "-1");

      string checkReplace;
      objectGet(editObj, "onlyReplaceCells", checkReplace, "all");
      if (checkReplace == "all")
      {
         edit.shouldReplaceCellTypes = false;
      }
      else
      {
         edit.shouldReplaceCellTypes = true;
         objectGet(editObj, "onlyReplaceCells", edit.replaceCellTypes);
      }
      
      string regionType;
      objectGet(editObj, "region", regionType, "");
      if (regionType == "ball")
      {
         vector<double> center;
         objectGet(editObj, "center", center);
         assert(center.size() == 3);
         double radius;
         objectGet(editObj, "radius", radius, "-1");
         assert(radius >= 0);
         edit.region = new Ball(center[0], center[1], center[2], radius);
      }
      else if (regionType == "box" || regionType=="ellipse")
      {
         vector<double> lower;
         objectGet(editObj, "lower", lower);
         assert(lower.size() == 3);
         vector<double> extent;
         objectGet(editObj, "extent", extent);
         assert(extent.size() == 3);
         for (int jj=0; jj<3; jj++)
         {
            assert(extent[jj] >= 0);
         }
         if (regionType == "box")
         {
            edit.region = new Box( lower[0],  lower[1],  lower[2],
                                  extent[0], extent[1], extent[2]);
         }
         else
         {
            edit.region = new EllipsoidInBox( lower[0],  lower[1],  lower[2],
                                             extent[0], extent[1], extent[2]);
         }
      }
      else
      {
         assert(0 && "Unrecognized region, I only accept box, ellipse, or radius");
      }
      edits.push_back(edit);
   }
   
   //open up the input bucket
   int nx, ny, nz;
   double dx, dy, dz;
   BucketOfBits* inBucket;
   OBJECT* hObj;
   PFILE* infile;
   {
      infile = Popen((indir+"/anatomy#").c_str(), "r", COMM_LOCAL);
      hObj = infile->headerObject;
      objectGet(hObj, "nx", nx, "0");
      objectGet(hObj, "ny", ny, "0");
      objectGet(hObj, "nz", nz, "0");
      objectGet(hObj, "dx", dx, "0");
      objectGet(hObj, "dy", dy, "0");
      objectGet(hObj, "dz", dz, "0");
      
      inBucket = readPioFile(infile);
   }
   
   int nRecords = inBucket->nRecords();
   int nFields = inBucket->nFields();
   int gidIndex = inBucket->getIndex("gid");
   int cellTypeIndex = inBucket->getIndex("cellType");
   int sigma11Index = inBucket->getIndex("sigma11");
   int sigma12Index = inBucket->getIndex("sigma12");
   int sigma13Index = inBucket->getIndex("sigma13");
   int sigma22Index = inBucket->getIndex("sigma22");
   int sigma23Index = inBucket->getIndex("sigma23");
   int sigma33Index = inBucket->getIndex("sigma33");
   assert(sigma11Index >= 0);
   assert(sigma12Index >= 0);
   assert(sigma13Index >= 0);
   assert(sigma22Index >= 0);
   assert(sigma23Index >= 0);
   assert(sigma33Index >= 0);

   const int lineSize = 10+5+6*12+1;
   
   PFILE* outfile = Popen((outdir+"/anatomy").c_str(), "w", COMM_LOCAL);
   if (myRank == 0)
   {
      string exeVersion;
      objectGet(hObj, "exe_version", exeVersion, "unknown");
      string creationDate;
      objectGet(hObj, "created_time", creationDate, "unknown");
      string endianKey;
      objectGet(hObj, "endian_key", endianKey, "875770417");

      time_t now = time(0);
      //write the header
      stringstream ss;
      ss <<
         "anatomy FILEHEADER {\n"
         "  exe_version = " << exeVersion << ";\n"
         "  created_time = " << creationDate << ";\n"
         "  datatype = FIXRECORDASCII;\n"
         "  nfiles = " << outfile->ngroup << ";\n"
         "  nrecord = " << nx*ny*nz <<";\n"
         "  lrec = " << lineSize << ";\n"
         "  endian_key = " << endianKey << ";\n"
         "  nfields = 8;\n"
         "  field_names = gid cellType sigma11 sigma12 sigma13 sigma22 sigma23 sigma33;\n"
         "  field_types = u u f f f f f f;\n"
         "  nx = " << nx << "; ny = " << ny << "; nz = " << nz << ";\n"
         "  field_units = 1 1 mS/mm mS/mm mS/mm mS/mm mS/mm mS/mm;\n"
         "  dx = " << dx << "; dy = " << dy << "; dz = " << dz << ";\n"
         "  modifiedAnatomy = "<< ctime(&now) << " ;\n"
         "}\n";

      string line = ss.str();
      Pwrite(line.c_str(), line.size(), 1, outfile);
   }
   
   //go through the gids and do the modification
   //for each gid
   for (int irecord=0; irecord<nRecords; irecord++)
   {
      BucketOfBits::Record inRecord = inBucket->getRecord(irecord);

      //get all the data so we can play with it.
      int gid;
      inRecord.getValue(gidIndex, gid);
      int cellType;
      inRecord.getValue(cellTypeIndex, cellType);
      double sigma11, sigma12, sigma13, sigma22, sigma23, sigma33;
      inRecord.getValue(sigma11Index, sigma11);
      inRecord.getValue(sigma12Index, sigma12);
      inRecord.getValue(sigma13Index, sigma13);
      inRecord.getValue(sigma22Index, sigma22);
      inRecord.getValue(sigma23Index, sigma23);
      inRecord.getValue(sigma33Index, sigma33);
      
      //get the xyz point for this gid
      int ix = gid % nx;
      int temp = gid / nx;
      int iy = temp % ny;
      int iz = temp / ny;
      const double xx = (ix+0.5)*dx;
      const double yy = (iy+0.5)*dy;
      const double zz = (iz+0.5)*dz;

      //for each region
      for (int iedit=0; iedit<edits.size(); iedit++)
      {
         //if the gid is within the region
         IRegion* region=edits[iedit].region;
         if (region->contains(xx,yy,zz))
         {
            if (edits[iedit].shouldReplaceCellTypes == true)
            {
               bool cellTypeInList = false;
               for (int icell=0; icell<edits[iedit].replaceCellTypes.size(); icell++)
               {
                  if (cellType == edits[iedit].replaceCellTypes[icell])
                  {
                     cellTypeInList = true;
                     break;
                  }
               }
               if (! cellTypeInList)
               {
                  continue;
               }
            }

            //do the cellType modifications
            if (edits[iedit].cellType >= 0)
            {
               cellType = edits[iedit].cellType;
            }

            //do the conductivity modifications
            if (edits[iedit].gil >= 0 && edits[iedit].git >= 0)
            {
               double gLongitudinal = edits[iedit].gil;
               double gTransverse = edits[iedit].git;
               double gNormal;
               if (edits[iedit].gin <= 0)
               {
                  gNormal = gTransverse;
               }
               else
               {
                  gNormal = edits[iedit].gin;
               }
               double sigmaT[3][3];
               sigmaT[0][0] = sigma11;
               sigmaT[0][1] = sigma12;
               sigmaT[0][2] = sigma13;
               sigmaT[1][0] = sigma12;
               sigmaT[1][1] = sigma22;
               sigmaT[1][2] = sigma23;
               sigmaT[2][0] = sigma13;
               sigmaT[2][1] = sigma23;
               sigmaT[2][2] = sigma33;

               //do the eigenvalue decomp
               double vT[3][3];
               double cond[3];
               symmetricTensorDecomp(sigmaT, vT, cond);
               double reltol = 1e-5;
               assert(cond[1]*(1+reltol) <= cond[2]);
               cond[2] = gLongitudinal;
               cond[1] = gTransverse;
               cond[0] = gNormal;

               //remake the tensor
               double result[3][3];
               for (int ii=0; ii<3; ii++) {
                  for (int kk=ii; kk<3; kk++) {
                     result[ii][kk] = 0;
                     for (int jj=0; jj<3; jj++) {
                        result[ii][kk] += vT[jj][ii]*cond[jj]*vT[jj][kk];
                     }
                  }
               }
               sigma11 = result[0][0];
               sigma12 = result[0][1];
               sigma13 = result[0][2];
               sigma22 = result[1][1];
               sigma23 = result[1][2];
               sigma33 = result[2][2];
            }
         }
      }
      //output data record.
      stringstream ss;
      ss << setw(10) << gid
         << setw(5) << cellType
         << " " << printToWidth(sigma11,11)
         << " " << printToWidth(sigma12,11)
         << " " << printToWidth(sigma13,11)
         << " " << printToWidth(sigma22,11)
         << " " << printToWidth(sigma23,11)
         << " " << printToWidth(sigma33,11)
         << endl
         ;
      
      string line = ss.str();
      assert(line.size()==lineSize);
      Pwrite(line.c_str(), line.size(), 1, outfile);
   }

   Pclose(infile);
   int retval = Pclose(outfile);
   return retval;
}
