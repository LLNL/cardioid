#include <map>
#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <iomanip>
#include <sstream>
#include <string>

#include <sys/stat.h>
#include <errno.h>

#include "Long64.hh"
#include "TupleToIndex.hh"

using namespace std;  // so sue me.

typedef map<pair<int, int>, int> DomainMap;



void readDomainInfo(istream& in, DomainMap& domainMap)
{
   string little, bunny, foo, hopping, thru, the, forest, bar;
   in>>little>>bunny>>foo>>foo>>hopping>>thru>>the>>forest>>foo>>foo>>bar ;

   int nDomains = 0;
   while (true)
   {
      int rank, nTissue, nBB, dz, dz4, nR, nD;
      double cost, timeR, timeD, timeRef;
      in >> rank >> nTissue >> nBB >> dz >> dz4 >> cost >> nR
         >> nD >> timeR >> timeD >> timeRef;
      if (in.eof()) break;
      ++nDomains;

      if (nTissue < 2) continue;
      pair<int, int> key = make_pair(nBB, nD);
      
      domainMap[key] = max(domainMap[key], nTissue);
   } 
   
   cout << "Read " << nDomains << " domains.  Kept "
        << domainMap.size()<< " unique representatives" << endl;
}

void writeAnatomy(ostream& out, int nx, int ny, int nz, int nTissue)
{
   int lRec = 21;
   assert (nTissue > 1);
   out <<
      "anatomy FILEHEADER {\n"
      "   datatype = FIXRECORDASCII;\n"
      "   nfiles = 1;\n"
      "   nfields = 4;\n"
      "   field_names = gid cellType theta phi;\n"
      "   field_types = u u u u;\n"
      "   nx = "<<nx+2<<"; ny = "<<ny+2<<"; nz = "<<nz+2<<";\n"
      "   nrecord = "<<nTissue<<";\n"
      "   lrec = "<<lRec<<";\n"
      "}\n\n";

   TupleToIndex tuple2index(nx+2, ny+2, nz+2);
   vector<Long64> allCells;
   allCells.reserve(nx*ny*nz);
   for (unsigned kk=1; kk<=nz; ++kk)
      for (unsigned jj=1; jj<=ny; ++jj)
         for (unsigned ii=1; ii<=nx; ++ii)
            allCells.push_back(tuple2index(ii, jj, kk));
                               
   out << setw(12) << allCells.back() << " 100 0 0\n";
   for (unsigned ii=0; ii<nTissue-1; ++ii)
      out <<setw(12) << allCells[ii] << " 100 0 0\n";
}

int DirTestCreate(const char *dirname)
{
   int mode = 0775;
   struct stat statbuf;
   int rc;
   rc = stat(dirname, &statbuf);
   if (rc == -1 && errno == ENOENT)
   {
      printf("Creating Directory: %s\n", dirname);
      rc = mkdir(dirname, mode);
      rc = stat(dirname, &statbuf);
   }
   if (rc != 0 || !(statbuf.st_mode & S_IFDIR))
   {
      printf("Can't Stat the Directory %s\n", dirname);
      printf("%d %x %x %x\n", rc, statbuf.st_mode, S_IFDIR, statbuf.st_mode & S_IFDIR);
   }
   return rc;
}

void writeInputDeck(ostream& out, int nD)
{
   out <<
      "simulate SIMULATE\n"
      "{\n"
      "   anatomy = pioAnatomy;\n"
      "   decomposition = koradi;\n"
      "   diffusion = FGR;\n"
      "   reaction = tt06dev;\n"
      "   stimulus = boxStimulus;\n"
      "   loop = 0;\n"            // in timesteps
      "   maxLoop = 10000;\n"    // in timesteps
      "   dt = 10 us;\n"
      "   time = 0;\n"            // msec
      "   printRate = 2500;\n"      // in timesteps
      "   snapshotRate = 1000001;\n" // in timesteps
      "   parallelDiffusionReaction = 1;\n"
      "   nDiffusionCores = "<<nD<<";\n"
      "}\n\n"
      "pioAnatomy ANATOMY\n"
      "{\n"
      "   method = pio;\n"
      "   fileName = snapshot.initial/anatomy#;\n"
      "   dx = 0.10;\n"   // in mm
      "   dy = 0.10;\n"   // in mm
      "   dz = 0.10;\n"   // in mm
      "   conductivity = conductivity;\n"
      "}\n\n"
      "koradi DECOMPOSITION\n"
      "{\n"
      "   method = koradi;\n"
      "   verbose = 1;\n"
      "   alpha = 0.25;\n"
      "   maxVoronoiSteps = 30;\n"
      "   maxSteps = 500;\n"
      "   tolerance = 0.01;\n"
      "   outputRate = 1000;\n"
      "   nbrDeltaR = 2;\n"
      "}\n\n"
      "FGR DIFFUSION\n"
      "{\n"
      "   method = FGR;\n"
      "   diffusionScale = 714.2857143;\n"      // mm^3/mF
      "}\n\n"
      "conductivity CONDUCTIVITY\n"
      "{\n"
      "    method = fibre;\n"
      "    sigmaLi = 0.0001334177;\n"   // units S/mm
      "    sigmaTi = 0.0000176062;\n"  // units S/mm
      "}\n\n"
      "tt06dev REACTION\n"
      "{\n"
      "    method = TT06Dev;\n"
      "    tolerance = 0.0001;\n"
      "    mod = 1;\n"
      "    fastGate =1;\n" 
      "    fastNonGate =1;\n" 
      "    cellTypes = endo mid epi;\n"
      "}\n\n"
      "endo CELLTYPE {clone=endoCellML;}\n"
      "mid  CELLTYPE {clone=midCellML;}\n"
      "epi  CELLTYPE {clone=epiCellML;}\n\n"
      "boxStimulus STIMULUS\n"
      "{\n"
      "   method = box;\n"
      "   vStim = -35.71429;\n"
      "   tStart = 0;\n"
      "   duration = 2;\n"
      "   period = 1000;\n"
      "}\n"
      ;
}



int main(int argc, char** argv)
{
   const int nx = 18;
   const int ny = 16;
   

   DomainMap domainMap;
   readDomainInfo(cin, domainMap);

   for (DomainMap::iterator
           iter = domainMap.begin(); iter!=domainMap.end(); ++iter)
   {
      int nBB = iter->first.first;
      int nD  = iter->first.second;
      int nTissue = iter->second;
      int nz = nBB/(nx*ny);

      stringstream dirname;
      dirname << nx<<"x"<<ny<<"x"<<nz<<"_"<<nTissue<<"/";
      string snapshotDir = dirname.str() +"/snapshot.initial/";
      DirTestCreate(dirname.str().c_str());
      DirTestCreate(snapshotDir.c_str());

      {
         string filename = snapshotDir +"anatomy#000000";
         ofstream aStream(filename.c_str());
         writeAnatomy(aStream, nx, ny, nz, nTissue);
      }
      {
         string filename = dirname.str() + "object.data";
         ofstream objstream(filename.c_str());
         writeInputDeck(objstream, nD);
      }
   }
   
}
