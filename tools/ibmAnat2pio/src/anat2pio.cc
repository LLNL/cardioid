// ToDo:

// 1. add data to header
//    - create time
//    - name of converted files
//
// 2.  Add binary file capability
//
// 3.  Proper parsing of command line with error checking
//
// 4.  Ability to convert only anatomy without orientation
//
// 5.  Filter out or convert cell types (such as type 0)
//
// 6.  Write cell types in header
//
// 7.  Add hmatrix
//
// 8.  Make sure all freads actually succeed.


#include <set>
#include <string>
#include <vector>
#include <cassert>
#include <iostream>

#include "EndianSwap.hh"
#include "heartIO.hh"
#include "TupleToIndex.hh"

using std::set;
using std::string;
using std::vector;
using std::endl;

#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE 1
#define _LARGE_FILE

typedef unsigned long long long64;



class CommandLineArguments
{
 public:
   CommandLineArguments(int argc, char** argv);

   string anatomyFile;
   string thetaFile;
   string phiFile;
   string asciiFile;
   string binaryFile;
};


int main(int argc, char** argv)
{
   CommandLineArguments options(argc, argv);

   FILE* anaFile    = fopen(options.anatomyFile.c_str(), "r");
   FILE* thetaFile  = fopen(options.thetaFile.c_str(), "r");
   FILE* phiFile    = fopen(options.phiFile.c_str(), "r");
   FILE* asciiFile  = fopen(options.asciiFile.c_str(), "w");
//   FILE* binaryFile = fopen(options.binaryFile.c_str(), "w");
   
   // get grid size info from anatomy file
   // Files as supplied from IBM are big endian.
   EndianSwap endianSwap(825373492);
   int nx, ny, nz;
   fseeko(anaFile, 26, SEEK_SET);
   fread(&nx, 4, 1, anaFile); endianSwap(nx);
   fread(&ny, 4, 1, anaFile); endianSwap(ny);
   fread(&nz, 4, 1, anaFile); endianSwap(nz);

   // ewd DEBUG
   //std::cout << "nx = " << nx << ", ny = " << ny << ", nz = " << nz << endl;

   
   string header = headSpace(2048);
   // leave space for header in pio files
   fwrite(header.c_str(), 1, header.size(), asciiFile);
//   fwrite(header.c_str(), 1, header.size(), binaryFile);

   vector<int> count(256, 0);
   
   // set each file stream to data position
   fseeko(anaFile,   274, SEEK_SET);
   fseeko(thetaFile, 274, SEEK_SET);
   fseeko(phiFile,   274, SEEK_SET);

   set<int> cellSet;
   long64 nRecA = 0;
   TupleToIndex tupleToIndex(nx, ny, nz);
   for (int kk=0; kk<nz; ++kk)
      for (int jj=0; jj<ny; ++jj)
	 for (int ii=0; ii<nx; ++ii)
	 {
	    unsigned char buf;
	    fread(&buf, 1, 1, anaFile);   int cellType = (int) buf;
	    fread(&buf, 1, 1, thetaFile); int theta    = (int) buf;
	    fread(&buf, 1, 1, phiFile);   int phi      = (int) buf;

	    cellSet.insert(cellType);
	    ++count[cellType];
	    Long64 gid = tupleToIndex(ii, jj, kk);

	    if (cellType != 0)
	    {
	       ++nRecA;
	       fprintf(asciiFile, "%11llu %3d %3d %3d\n",
		       gid, cellType, theta, phi);
	    }
	    
	 }

   fclose(anaFile);
   fclose(thetaFile);
   fclose(phiFile);
   fclose(asciiFile);
//   fclose(binaryFile);
   
   HeaderInfo aHdr;
   aHdr.dataType = "FIXRECORDASCII";
   aHdr.nFiles = 1;
   aHdr.nRec = nRecA;
   aHdr.lRec = 24;
   aHdr.nFields = 4;
   aHdr.fieldNames = "gid cellType theta phi";
   aHdr.fieldTypes = "u u u u";
   HeaderInfo bHdr = aHdr;
   bHdr.dataType = "FIXRECORDBINARY";
   
   writeHeader(options.asciiFile,  aHdr, nx, ny, nz, cellSet);
//   writeHeader(options.binaryFile, bHdr, nx, ny, nz, cellSet);


   for (set<int>::iterator iter=cellSet.begin(); iter!=cellSet.end(); ++iter)
      printf("%d: %d\n", *iter, count[*iter]);
   
}

CommandLineArguments::CommandLineArguments(int argc, char** argv)
{
   anatomyFile = argv[1];
   thetaFile   = argv[2];
   phiFile     = argv[3];

   asciiFile = "anatomyAscii#000000";
   binaryFile = "anatomyBinary#000000";
}
   
