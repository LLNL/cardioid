#include "Array3d.hh"

int main(int argc, char** argv)
{
   Array3d<int> aa;
   int nx = 3;
   int ny = 5;
   int nz = 7;
   aa.resize(nx, ny, nz);


   int data = 32;

   for (int jj=0; jj<ny; ++jj)
      for (int ii=0; ii<nx; ++ii)
	 for (int kk=0; kk<nz; ++kk)
	 {
	    aa(ii, jj, kk) = data;
	    assert( aa(ii, jj, kk) == data);
	    assert( aa.cArray()[ii][jj][kk] == data);
	    data++;
	 }

   for (unsigned ii=0; ii<aa.size(); ++ii)
      aa(ii) = 0;

   for (int jj=0; jj<ny; ++jj)
      for (int ii=0; ii<nx; ++ii)
	 for (int kk=0; kk<nz; ++kk)
	    assert(aa.cArray()[ii][jj][kk] == 0);

   for (int jj=0; jj<ny; ++jj)
      for (int ii=0; ii<nx; ++ii)
	 for (int kk=0; kk<nz; ++kk)
	 {
	    
	    aa.cArray()[ii][jj][kk] = data;
	    assert( aa(ii, jj, kk) == data);
	    assert( aa(ii, jj, kk)  == data);
	    data++;
	 }
}

   
