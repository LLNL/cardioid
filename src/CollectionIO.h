#ifndef COLLECTIONIO_H
#define COLLECTIONIO_H

#include <iostream>
#include <iomanip>
using namespace std;

class CollectionIO
{

  public:
  
  CollectionIO();
  ~CollectionIO();
  
  int readHeader();
  int readDataAscii(int nlines);
  int readDataBinary(int nrec);
  
};
#endif

