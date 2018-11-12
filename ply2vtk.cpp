// takes in vtk and generates list of x,y,z locations at centroid of each tetrahedral element

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

using namespace std;

int main(int argc, char **argv) 
{
// READ IN PLY
   
   string dummy;
   
   ifstream fin;
   int numpoints, numcells;
   
   fin.open(argv[1]);
   getline(fin, dummy);
   getline(fin, dummy);
   getline(fin, dummy);
   getline(fin, dummy);
   fin >> dummy >> dummy >> numpoints >> dummy;
   getline(fin, dummy);
   getline(fin, dummy);
   getline(fin, dummy);   
   fin >> dummy >> dummy >> numcells >> dummy;  
   getline(fin, dummy);
   getline(fin, dummy);   
   
   int no_of_cols = 3;
   int no_of_rows = numpoints;
   int initial_value = 0;

   vector<vector<double> > points;
   points.resize(no_of_rows, vector<double>(no_of_cols, initial_value));
   
   for (int i = 0; i < numpoints; i++){
     fin >> points[i][0] >> points[i][1] >> points[i][2];
   }
   
   no_of_cols = 4;
   no_of_rows = numcells;
   initial_value = 0;

   vector<vector<int> > triverts;
   triverts.resize(no_of_rows, vector<int>(no_of_cols, initial_value));
   
   for (int i = 0; i < numcells; i++){
      fin >> triverts[i][0] >> triverts[i][1] >> triverts[i][2] >> triverts[i][3];
   }
   
   fin.close();
   
   
// WRITE OUT VTK UNSTRUCTURED GRID FILE
   
   ofstream fout;
   fout.open(argv[2]); 
   fout << "# vtk DataFile Version 3.0" << endl;
   fout << "vtk output" << endl;
   fout << "ASCII" << endl;
   fout << "DATASET UNSTRUCTURED_GRID" << endl;
   fout << "POINTS " << numpoints << " double" << endl;
      for (int i = 0; i < numpoints; i++){
      fout << points[i][0] << " " << points[i][1] << " " << points[i][2] << endl;
   }
   
   fout << "CELLS " << numcells << " " << 4 * numcells << endl;
   for (int i = 0; i < numcells; i++){
      fout << triverts[i][0] << " " << triverts[i][1] << " " << triverts[i][2] << " " << triverts[i][3] << endl;
   }
   
   fout << "CELL_TYPES " << numcells << endl;
   for (int i = 0; i < numcells; i++){
      fout << "5" << endl;
   }
   
   fout.close();
   
   return 0;
}
