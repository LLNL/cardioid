
#include "heartIO.hh"
#include "Long64.hh"
#include "TupleToIndex.hh"

#include <string>
#include <set>
#include <cmath>
#include <iostream>

using namespace std;


int xyTorus(FILE* file, int nx, int ny, int nz);
int solidBrick(FILE* file, int nx, int ny, int nz);
int paraboloid(FILE* file, int nx, int ny, int nz);




int main(int argc, char** argv)
{

   int nx=100;
   int ny=100;
   int nz=3;


   string filename = "anatomy#000000";

   FILE* file = fopen(filename.c_str(), "w");
   
   string header = headSpace(2048);
   fwrite(header.c_str(), 1, header.size(), file);

//   int nRecords = xyTorus(file, nx, ny, nz);
//    nx=ny=nz=310;
//    int nRecords = solidBrick(file, nx, ny, nz);
   nx=ny=400; nz=800;
   int nRecords = paraboloid(file, nx, ny, nz);

   fclose(file);

   HeaderInfo hdr;
   hdr.dataType = "FIXRECORDASCII";
   hdr.nFiles = 1;
   hdr.nRec = nRecords;
   hdr.lRec = 24;
   hdr.fieldNames = "gid cellType theta phi";
   hdr.fieldTypes = "u u u u";
      
   
   writeHeader(filename, hdr, nx, ny, nz, set<int>());

}

int xyTorus(FILE* file, int nx, int ny, int nz)
{
   double xCenter = double(nx-1)/2.0;
   double yCenter = double(ny-1)/2.0;
   double zCenter = double(nz-1)/2.0;

   TupleToIndex tupleToIndex(nx, ny, nz);
   int nRecords = 0;
   for (int kk=0; kk<nz; ++kk)
   {
      double rz = kk - zCenter;
      for (int jj=0; jj<ny; ++jj)
      {
	 double ry = jj - yCenter;
	 for (int ii=0; ii<nx; ++ii)
	 {
	    double rx = ii-xCenter;
	    double r = sqrt(rx*rx + ry*ry + rz*rz);

	    int cellType = 100;
	    if (r>50)
	       cellType = 0;
	    if (r<20)
	       cellType = 9;

	    ++nRecords;
	    Long64 gid = tupleToIndex(ii, jj, kk);
	    int theta = 0;
	    int phi = 0;

	    fprintf(file, "%11llu %3d %3d %3d\n", gid, cellType, theta, phi);
	 }
      }
   }
   return nRecords;
}

int solidBrick(FILE* file, int nx, int ny, int nz)
{
   double xCenter = double(nx-1)/2.0;
   double yCenter = double(ny-1)/2.0;
   double zCenter = double(nz-1)/2.0;

   TupleToIndex tupleToIndex(nx, ny, nz);
   int nRecords = 0;
   for (int kk=0; kk<nz; ++kk)
   {
      double rz = kk - zCenter;
      for (int jj=0; jj<ny; ++jj)
      {
	 double ry = jj - yCenter;
	 for (int ii=0; ii<nx; ++ii)
	 {
	    double rx = ii-xCenter;

	    int cellType = 100;

	    ++nRecords;
	    Long64 gid = tupleToIndex(ii, jj, kk);
	    int theta = 0;
	    int phi = 0;

	    fprintf(file, "%11llu %3d %3d %3d\n", gid, cellType, theta, phi);
	 }
      }
   }
   return nRecords;
}

double paraboloidRadius(int rMax, int nz, int z)
{
   if (z==0)
      return 0;
   if (z>nz/2)
      return rMax;

   double a = nz/(2.*rMax*rMax);
   
   return sqrt(z/a);   
}


// This isn't actually a paraboloid, just a rough approximation.
int paraboloid(FILE* file, int nx, int ny, int nz)
{
   double xCenter = double(nx-1)/2.0;
   double yCenter = double(ny-1)/2.0;
   double zCenter = double(nz-1)/2.0;

   TupleToIndex tupleToIndex(nx, ny, nz);
   int nRecords = 0;
   for (int kk=0; kk<nz; ++kk)
   {
      double rz = kk - zCenter;
      double rMax = paraboloidRadius(nx/2, nz, kk);
      double rMin = max(0., rMax - 130);
      for (int jj=0; jj<ny; ++jj)
      {
	 double ry = jj - yCenter;
	 for (int ii=0; ii<nx; ++ii)
	 {
	    double rx = ii-xCenter;

	    double r = sqrt(rx*rx + ry*ry );

	    int cellType = 100;
	    if (r>rMax)
	       cellType = 0;
	    if (r<rMin)
	       cellType = 9;
	    
	    ++nRecords;
	    Long64 gid = tupleToIndex(ii, jj, kk);
	    int theta = 0;
	    int phi = 0;

	    fprintf(file, "%11llu %3d %3d %3d\n", gid, cellType, theta, phi);
	 }
      }
   }
   return nRecords;
   
}

