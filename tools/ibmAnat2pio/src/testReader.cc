#include <string>
#include <vector>
#include <cassert>
#include <iostream>
using namespace std;

#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE 1
#define _LARGE_FILE

typedef unsigned long long long64;

int main(int argc, char** argv)
{

  if (argc != 2)
  {
    cout << "Usage:  testReader [anatomy filename]" << endl;
    exit(1);
  }
  string anatomyFile = argv[1];
  FILE* anaFile    = fopen(anatomyFile.c_str(), "r");
   
  // get grid size info from anatomy file
  // Files as supplied are big endian.
  int nx, ny, nz;
  fseeko(anaFile, 26, SEEK_SET);
  fread(&nx, 4, 1, anaFile);
  fread(&ny, 4, 1, anaFile);
  fread(&nz, 4, 1, anaFile);

  // ewd DEBUG
  cout << "nx = " << nx << ", ny = " << ny << ", nz = " << nz << endl;

  fseeko(anaFile, 274, SEEK_SET);

  int size = nx*ny*nz;
  int data;
  for (int i=0; i<size; i++)
  {
    fread(&data, sizeof(int), 1, anaFile);
    cout << "i = " << i << ", data = " << data << endl;
  }
  

  fclose(anaFile);
   
}

