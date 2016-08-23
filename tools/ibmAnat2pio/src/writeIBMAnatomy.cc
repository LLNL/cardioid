#include <string>
#include <vector>
#include <cassert>
#include <iostream>
using namespace std;

#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE 1
#define _LARGE_FILE

typedef unsigned long long long64;


class CommandLineArguments
{
 public:
    CommandLineArguments(int argc, char** argv);
    int nx, ny, nz;
    string binaryFile;
    string thetaFile;
    string phiFile;
};


int main(int argc, char** argv)
{
  int celltype = 100;  // 100 = TT04_Endo, 101 = TT04_Mid, 102 = TT04_Epi

  CommandLineArguments options(argc, argv);
  FILE* binaryFile = fopen(options.binaryFile.c_str(), "w");
  FILE* thetaFile = fopen(options.thetaFile.c_str(), "w");
  FILE* phiFile = fopen(options.phiFile.c_str(), "w");
   
  int grid[3];
  grid[0] = options.nx;
  grid[1] = options.ny;
  grid[2] = options.nz;

  // write header:  22 bytes, four ints, 274 byte header total
  char tmp[1] = {'x'};
  for (int i=0; i<22; i++)
  {
    fwrite(&tmp[0], 1, sizeof(tmp), binaryFile);
    fwrite(&tmp[0], 1, sizeof(tmp), thetaFile);
    fwrite(&tmp[0], 1, sizeof(tmp), phiFile);
  }
  
  // temporary int
  fwrite(&grid[0], 1, sizeof(int), binaryFile);
  fwrite(&grid[0], 1, sizeof(int), thetaFile);
  fwrite(&grid[0], 1, sizeof(int), phiFile);

  // write grid size
  fwrite(&grid[0], 3, sizeof(int), binaryFile);
  fwrite(&grid[0], 3, sizeof(int), thetaFile);
  fwrite(&grid[0], 3, sizeof(int), phiFile);

  // pad out header to 274 bytes total
  for (int i=38; i<274; i++)
  {
    fwrite(&tmp[0], 1, sizeof(tmp), binaryFile);
    fwrite(&tmp[0], 1, sizeof(tmp), thetaFile);
    fwrite(&tmp[0], 1, sizeof(tmp), phiFile);
  }
  
  // fill anatomy data array uniformly, with zeroes at edges
  int size = grid[0]*grid[1]*grid[2];
  int anatdata[size];
  int cnt = 0;
  for (int i=0; i<grid[0]; i++)
    for (int j=0; j<grid[1]; j++)
      for (int k=0; k<grid[2]; k++)
      {
        int val = celltype;
        if (i == 0 || i == grid[0]-1 || j == 0 || j == grid[1]-1 || k == 0 || k == grid[2]-1)
          val = 0;
        anatdata[cnt++] = val;
      }
  // write anatomy data
  fwrite(&anatdata[0], size, sizeof(int), binaryFile);

  // fill theta, phi data arrays
  int thetadata[size];
  int phidata[size];
  for (int i=0; i<size; i++)
  {
    thetadata[i] = i%358;
    phidata[i] = i%358;
  }
  
  // write angle data
  fwrite(&thetadata[0], size, sizeof(int), thetaFile);
  fwrite(&phidata[0], size, sizeof(int), phiFile);
  
  
  fclose(binaryFile);
  fclose(thetaFile);
  fclose(phiFile);
   
}

CommandLineArguments::CommandLineArguments(int argc, char** argv)
{
  if (argc != 7)
  {
    cout << "Usage:  writeIBMAnatomy [nx] [ny] [nz] [anatomy filename] [theta filename] [phi filename]" << endl;
    exit(1);
  }
  nx = atoi(argv[1]);
  ny = atoi(argv[2]);
  nz = atoi(argv[3]);
  binaryFile = argv[4];
  thetaFile = argv[5];
  phiFile = argv[6];
}
   
