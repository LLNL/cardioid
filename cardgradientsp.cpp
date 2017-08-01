#include "mfem.hpp"
#include "cardfiber.h"
#include "genfiber.h"
#include "constants.h"
#include "triplet.h"
#include "io.h"
#include "kdtree++/kdtree.hpp"
#include "cardgradientsp.h"
#include <sstream>
#include <iomanip>
#include "pio.h"
#include "ioUtils.h"
#include "heap.h"

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

using namespace std;
using namespace mfem;

void getCardGradientsp(Mesh* mesh, GridFunction& x_psi_ab, GridFunction& x_phi_epi, GridFunction& x_phi_lv, GridFunction& x_phi_rv,
        tree_type& kdtree, vector<vector<int> >& vert2Elements, vector<Vector>& boundingbox, double dd,
        Vector& conduct, Vector& fiberAngles, int num_procs, int myid) {
    Vector min = boundingbox[0];
    Vector max = boundingbox[1];

    double xmin = min(0);
    double ymin = min(1);
    double zmin = min(2);

    double x_dim = max(0) - xmin;
    double y_dim = max(1) - ymin;
    double z_dim = max(2) - zmin;

    double dd10 = dd * 10;

    const double dx = x_dim / (int(x_dim / (dd10))*10 - 1);
    const double dy = y_dim / (int(y_dim / (dd10))*10 - 1);
    const double dz = z_dim / (int(z_dim / (dd10))*10 - 1);

    int nx = int(x_dim / dx) + 1;
    int ny = int(y_dim / dy) + 1;
    int nz = int(z_dim / dz) + 1;

    int totalCardPoints = 0;
    vector<anatomy> anatVectors;
    
    double cutoff=getMaxEdgeLen(mesh)*0.6123724356957945;  //Radius of circumsphere sqrt(6)/4 
    if(myid==0){
        cout << "\nCutoff for nearest point is " << cutoff << std::endl;    
    }
    
    long long gid_dim = nx * ny*nz;
    // MPI Parallel
    for (long long g = myid; g < gid_dim; g += num_procs) {
        int i = g % nx;
        int j = (g / nx) % ny;
        int k = g / nx / ny;

        double x = xmin + i*dx;
        double y = ymin + j*dy;
        double z = zmin + k*dz;
        triplet pt(x, y, z, 0);
        std::pair<tree_type::const_iterator, double> found = kdtree.find_nearest(pt);
        assert(found.first != kdtree.end());
        // Skip if the distance between pt and nearest is larger than cutoff
        if (found.second>cutoff) continue; 
         
        //For barycentric
        Vector q(4);
        q(0) = x;
        q(1) = y;
        q(2) = z;
        q(3) = 1.0;
        
        triplet vetexNearPt = *found.first;
        int vertex = vetexNearPt.getIndex();
        
        vector<int> elements = vert2Elements[vertex];
        for (unsigned e = 0; e < elements.size(); e++) {
            int eleIndex = elements[e];

            if (isInTetElement(q, mesh, eleIndex)) {                
                anatomy anat;
                anat.gid = g;
                calcGradient(x_psi_ab, x_phi_epi, x_phi_lv, x_phi_rv, conduct, fiberAngles, q, eleIndex, anat);                
                anatVectors.push_back(anat);

                totalCardPoints++;
                if (totalCardPoints % 10000 == 0) {
                    cout << "Processor " << myid << " finish " << totalCardPoints << " points." << endl;
                    cout.flush();
                }
                break; // If the point is found in an element, don't need to check next one in the list. 
            }
        }

    }
    filerheader header;
    header.nx = nx;
    header.ny = ny;
    header.nz = nz;
    header.dx = dx;
    header.dy = dy;
    header.dz = dz;

    int globalTotCardPoints;

    MPI_Allreduce(&totalCardPoints, &globalTotCardPoints, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    header.nrecord = globalTotCardPoints;

    // Parallel I/O
    string fullname = "snapshot.initial";
    if (myid == 0) {
        DirTestCreate(fullname.c_str());
    }
    
    fullname += "/anatomy";
    int lrec = 80;
    heap_allocate(lrec*totalCardPoints*64 + 4096);
    
    PFILE* file = Popen(fullname.c_str(), "w", MPI_COMM_WORLD);
    PioReserve(file, lrec*totalCardPoints*64 + 4096);
    //int nFiles_ = num_procs;
//    if (nFiles_ > 0)
//        PioSet(file, "ngroup", nFiles_);  

    if (myid == 0) {
        Pprintf(file, "anatomy FILEHEADER { \n");
        Pprintf(file, "  datatype = VARRECORDASCII;\n");
        Pprintf(file, "  nfiles = %d;  \n", file.nfiles);
        Pprintf(file, "  nrecord = %d; \n", header.nrecord);
        Pprintf(file, "  nfields = 8; \n");
        Pprintf(file, "  lrec = %d; \n", lrec);
        Pprintf(file, "  endian_key = 875770417; \n");
        Pprintf(file, "  field_names = gid cellType sigma11 sigma12 sigma13 sigma22 sigma23 sigma33; \n");
        Pprintf(file, "  field_types = u u f f f f f f; \n");
        Pprintf(file, "  field_units = 1 1 mS/mm mS/mm mS/mm mS/mm mS/mm mS/mm; \n");
        Pprintf(file, "  nx =  %d;\n", header.nx);
        Pprintf(file, "  ny =  %d;\n", header.ny);
        Pprintf(file, "  nz =  %d;\n", header.nz);
        Pprintf(file, "  dx =  %f;\n", header.dx);
        Pprintf(file, "  dy =  %f;\n", header.dy);
        Pprintf(file, "  dz =  %f;\n", header.dz);
        Pprintf(file, "} \n\n");

    }  

    for (unsigned i = 0; i < anatVectors.size(); i++) {
        anatomy anat = anatVectors[i];
        Pprintf(file, "    %llu  %d ", anat.gid, anat.celltype);
        for (int j = 0; j < 6; j++) {
            Pprintf(file, "%f ", anat.sigma[j]);
        }
        Pprintf(file, "\n");
    }

    Pclose(file);
    
}

template < class ContainerT >
void tokenize(const std::string& str, ContainerT& tokens,
              const std::string& delimiters = " ", const bool trimEmpty = true) {
    typedef ContainerT Base;
    typedef typename Base::value_type ValueType;
    typedef typename ValueType::size_type SizeType;
    SizeType pos, lastPos = 0;
    while (true) {
        pos = str.find_first_of(delimiters, lastPos);
        if (pos == std::string::npos) {
            pos = str.length();

            if (pos != lastPos || !trimEmpty)
                tokens.push_back(ValueType(str.data() + lastPos, (SizeType)pos - lastPos));

            break;
        } else {
            if (pos != lastPos || !trimEmpty)
                tokens.push_back(ValueType(str.data() + lastPos, (SizeType)pos - lastPos));
        }

        lastPos = pos + 1;
    }
};

void readlines(MPI_File *in, const int rank, const int size, const int overlap,
               char ***lines, int *nlines) {
    MPI_Offset filesize;
    MPI_Offset localsize;
    MPI_Offset start;
    MPI_Offset end;
    char *chunk;

    /* figure out who reads what */

    MPI_File_get_size(*in, &filesize);
    localsize = filesize/size;
    start = rank * localsize;
    end   = start + localsize - 1;

    //printf("1. Process %d start= %d end= %d\n", rank, start, end);

    /* add overlap to the end of everyone's chunk... */
    end += overlap;

    /* except the last processor, of course */
    if (rank == size-1) end = filesize;

    localsize =  end - start + 1;
    //printf("2. Process %d start= %d end= %d localsize=%d\n", rank, start, end, localsize);

    /* allocate memory */
    chunk = (char*) malloc( (localsize + 1)*sizeof(char));

    /* everyone reads in their part */
    MPI_File_read_at_all(*in, start, chunk, localsize, MPI_CHAR, MPI_STATUS_IGNORE);
    chunk[localsize] = '\0';

    /*
     * everyone calculate what their start and end *really* are by going 
     * from the first newline after start to the first newline after the
     * overlap region starts (eg, after end - overlap + 1)
     */

    int locstart=0, locend=localsize;
    if (rank != 0) {
        while(chunk[locstart] != '\n') locstart++;
        locstart++;
    }
    if (rank != size-1) {
        locend-=overlap;
        while(chunk[locend] != '\n') locend++;
    }
    localsize = locend-locstart+1;

    /* Now let's copy our actual data over into a new array, with no overlaps */
    char *data = (char *)malloc((localsize+1)*sizeof(char));
    memcpy(data, &(chunk[locstart]), localsize);
    data[localsize] = '\0';
    free(chunk);

    /* Now we'll count the number of lines */
    *nlines = 0;
    for (int i=0; i<localsize; i++)
        if (data[i] == '\n') (*nlines)++;

    /* Now the array lines will point into the data array at the start of each line */
    /* assuming nlines > 1 */
    *lines = (char **)malloc((*nlines)*sizeof(char *));
    (*lines)[0] = strtok(data,"\n");
    for (int i=1; i<(*nlines); i++)
        (*lines)[i] = strtok(NULL, "\n");

    return;
}

void getRotMatrixp(Mesh* mesh, GridFunction& x_psi_ab, GridFunction& x_phi_epi, GridFunction& x_phi_lv, GridFunction& x_phi_rv,
                   tree_type& kdtree, vector<vector<int> >& vert2Elements, Vector& fiberAngles, const char *fiblocs, int size, int rank)
{

   long long totalCardPoints = 0;
   
   //f_ofs << "# elementnum mat11 mat12 mat13 mat21 mat22 mat23 mat31 mat32 mat33" << endl;


   MPI_File in;

   int ierr = MPI_File_open(MPI_COMM_WORLD, (char *)fiblocs, MPI_MODE_RDONLY, MPI_INFO_NULL, &in);
   if (ierr)
   {
      if (rank == 0) fprintf(stderr, "Couldn't open file %s\n", fiblocs);
      MPI_Finalize();
      exit(2);
   }

   const int overlap = 200;
   char **lines;
   int nlines;
   readlines(&in, rank, size, overlap, &lines, &nlines);
   printf("Rank %d has %d lines\n", rank, nlines);

   std::string fileLine;

   const std::string comment = "#";

   vector<string> outLines;

   for (int i = 0; i < nlines; i++)
   {
      fileLine = lines[i];
      if (fileLine.compare(0, 1, comment) == 0) continue;
      std::vector<std::string> tokens;
      tokenize(fileLine, tokens);
      if (tokens.size() > 3)
      {
         double x = atof(tokens[1].c_str());
         double y = atof(tokens[2].c_str());
         double z = atof(tokens[3].c_str());
         triplet pt(x, y, z, 0);
         std::pair<tree_type::const_iterator, double> found = kdtree.find_nearest(pt);
         assert(found.first != kdtree.end());
         triplet vetexNearPt = *found.first;
         int vertex = vetexNearPt.getIndex();
         vector<int> elements = vert2Elements[vertex];
         bool findPt=false;
         for (unsigned e = 0; e < elements.size(); e++)
         {
            int eleIndex = elements[e];
            //For barycentric
            Vector q(4);
            q(0) = x;
            q(1) = y;
            q(2) = z;
            q(3) = 1.0;
            vector<double> barycentric;
            if (isInTetElement(q, mesh, eleIndex))
            {
               //cout << "fiblocs element index=" << tokens[0] << "; k-D tree index=" << eleIndex << endl;
               Vector psi_ab_vec(3);
               double psi_ab = 0.0;
               getCardEleGrads(x_psi_ab, q, eleIndex, psi_ab_vec, psi_ab);

               Vector phi_epi_vec(3);
               double phi_epi = 0.0;
               getCardEleGrads(x_phi_epi, q, eleIndex, phi_epi_vec, phi_epi);

               Vector phi_lv_vec(3);
               double phi_lv = 0.0;
               getCardEleGrads(x_phi_lv, q, eleIndex, phi_lv_vec, phi_lv);

               Vector phi_rv_vec(3);
               double phi_rv = 0.0;
               getCardEleGrads(x_phi_rv, q, eleIndex, phi_rv_vec, phi_rv);

               DenseMatrix QPfib(dim3, dim3);
               biSlerpCombo(QPfib, psi_ab, psi_ab_vec, phi_epi, phi_epi_vec,
                            phi_lv, phi_lv_vec, phi_rv, phi_rv_vec, fiberAngles);

               stringstream f_ofs;
               f_ofs << tokens[0] << " ";
               for (int ii = 0; ii < dim3; ii++)
               {
                  for (int jj = 0; jj < dim3; jj++)
                  {
                     f_ofs << QPfib(ii, jj) << " ";
                  }
               }
               f_ofs << endl;
               outLines.push_back(f_ofs.str());

               totalCardPoints++;
                if (totalCardPoints % 10000 == 0) {
                    cout << "Processor " << rank << " finish " << totalCardPoints << " points." << endl;
                    cout.flush();
                }
               findPt=true;
               break; // If the point is found in an element, don't need to check next one in the list. 

            }
         }
         
         if(!findPt){
               stringstream f_ofs;
               f_ofs << tokens[0] << " ";
               for (int ii = 0; ii < dim3; ii++)
               {
                  for (int jj = 0; jj < dim3; jj++)
                  {
                     f_ofs << "999 ";
                  }
               }
               f_ofs << endl;
               outLines.push_back(f_ofs.str());            
         }
         
         
      }
   }
   
   cout << "Processor " << rank << " has " << outLines.size() << " lines." << endl;

//    // Parallel I/O
//    string fullname = "omar";
//    if (rank == 0) {
//        DirTestCreate(fullname.c_str());
//    }
//    
//    fullname += "/rotmatrix";
//    int lrec = 80;
//    heap_allocate(lrec*totalCardPoints*64 + 4096);
//    
//    PFILE* file = Popen(fullname.c_str(), "w", MPI_COMM_WORLD);
//    PioReserve(file, lrec*totalCardPoints*64 + 4096);
//
//    for (unsigned i = 0; i < outLines.size(); i++) {
//        string line = outLines[i];
//        Pprintf(file, "%s", line.c_str());
//    }
//
//    Pclose(file);   
   
   int file_free = 0;
   MPI_Status status;

   if (rank == 0)
   {
      file_free = 1;
   }
   else
   {
      MPI_Recv(&file_free, 1, MPI_INT, rank - 1, 1, MPI_COMM_WORLD, &status);
   }

   if (file_free == 1)
   {
      ofstream out;
      
      if (rank == 0)
      {
         out.open("rotmatrix.txt");
         out << "# elementnum mat11 mat12 mat13 mat21 mat22 mat23 mat31 mat32 mat33" << endl;

      }
      else
      {
         out.open("rotmatrix.txt", std::fstream::app);
      }
      
      for (int ii = 0; ii < outLines.size(); ii++)
      {
         out << outLines[ii];
      }
      out.close();

   }

   if (rank != size - 1)
   {
      MPI_Send(&file_free, 1, MPI_INT, rank + 1, 1, MPI_COMM_WORLD);
   }


}

void getRotMatrixFastp(Mesh* mesh, GridFunction& x_psi_ab, GridFunction& x_phi_epi, GridFunction& x_phi_lv, GridFunction& x_phi_rv,
                   vector<vector<int> >& vert2Elements, Vector& fiberAngles, const char *fiblocs, int size, int rank)
{

   long long totalCardPoints = 0;
   
   //f_ofs << "# elementnum mat11 mat12 mat13 mat21 mat22 mat23 mat31 mat32 mat33" << endl;


   MPI_File in;

   int ierr = MPI_File_open(MPI_COMM_WORLD, (char *)fiblocs, MPI_MODE_RDONLY, MPI_INFO_NULL, &in);
   if (ierr)
   {
      if (rank == 0) fprintf(stderr, "Couldn't open file %s\n", fiblocs);
      MPI_Finalize();
      exit(2);
   }

   const int overlap = 200;
   char **lines;
   int nlines;
   readlines(&in, rank, size, overlap, &lines, &nlines);
   printf("Rank %d has %d lines\n", rank, nlines);

   std::string fileLine;

   const std::string comment = "#";

   vector<string> outLines;

   for (int i = 0; i < nlines; i++)
   {
      fileLine = lines[i];
      if (fileLine.compare(0, 1, comment) == 0) continue;
      std::vector<std::string> tokens;
      tokenize(fileLine, tokens);
      if (tokens.size() > 3)
      {
         int eleIndex=atoi(tokens[0].c_str());
         double x = atof(tokens[1].c_str());
         double y = atof(tokens[2].c_str());
         double z = atof(tokens[3].c_str());

            //For barycentric
            Vector q(4);
            q(0) = x;
            q(1) = y;
            q(2) = z;
            q(3) = 1.0;
            vector<double> barycentric;
            if (isInTetElement(q, mesh, eleIndex))
            {
               //cout << "fiblocs element index=" << tokens[0] << "; k-D tree index=" << eleIndex << endl;
               Vector psi_ab_vec(3);
               double psi_ab = 0.0;
               getCardEleGrads(x_psi_ab, q, eleIndex, psi_ab_vec, psi_ab);

               Vector phi_epi_vec(3);
               double phi_epi = 0.0;
               getCardEleGrads(x_phi_epi, q, eleIndex, phi_epi_vec, phi_epi);

               Vector phi_lv_vec(3);
               double phi_lv = 0.0;
               getCardEleGrads(x_phi_lv, q, eleIndex, phi_lv_vec, phi_lv);

               Vector phi_rv_vec(3);
               double phi_rv = 0.0;
               getCardEleGrads(x_phi_rv, q, eleIndex, phi_rv_vec, phi_rv);

               DenseMatrix QPfib(dim3, dim3);
               biSlerpCombo(QPfib, psi_ab, psi_ab_vec, phi_epi, phi_epi_vec,
                            phi_lv, phi_lv_vec, phi_rv, phi_rv_vec, fiberAngles);

               stringstream f_ofs;
               f_ofs << tokens[0] << " ";
               for (int ii = 0; ii < dim3; ii++)
               {
                  for (int jj = 0; jj < dim3; jj++)
                  {
                     f_ofs << QPfib(ii, jj) << " ";
                  }
               }
               f_ofs << endl;
               outLines.push_back(f_ofs.str());

               totalCardPoints++;
                if (totalCardPoints % 10000 == 0) {
                    cout << "Processor " << rank << " finish " << totalCardPoints << " points." << endl;
                    cout.flush();
                }


            
         }
         

         
         
      }
   }
   
   cout << "Processor " << rank << " has " << outLines.size() << " lines." << endl;

//    // Parallel I/O
//    string fullname = "omar";
//    if (rank == 0) {
//        DirTestCreate(fullname.c_str());
//    }
//    
//    fullname += "/rotmatrix";
//    int lrec = 80;
//    heap_allocate(lrec*totalCardPoints*64 + 4096);
//    
//    PFILE* file = Popen(fullname.c_str(), "w", MPI_COMM_WORLD);
//    PioReserve(file, lrec*totalCardPoints*64 + 4096);
//
//    for (unsigned i = 0; i < outLines.size(); i++) {
//        string line = outLines[i];
//        Pprintf(file, "%s", line.c_str());
//    }
//
//    Pclose(file);   
   
   int file_free = 0;
   MPI_Status status;

   if (rank == 0)
   {
      file_free = 1;
   }
   else
   {
      MPI_Recv(&file_free, 1, MPI_INT, rank - 1, 1, MPI_COMM_WORLD, &status);
   }

   if (file_free == 1)
   {
      ofstream out;
      
      if (rank == 0)
      {
         out.open("rotmatrix.txt");
         out << "# elementnum mat11 mat12 mat13 mat21 mat22 mat23 mat31 mat32 mat33" << endl;

      }
      else
      {
         out.open("rotmatrix.txt", std::fstream::app);
      }
      
      for (int ii = 0; ii < outLines.size(); ii++)
      {
         out << outLines[ii];
      }
      out.close();

   }

   if (rank != size - 1)
   {
      MPI_Send(&file_free, 1, MPI_INT, rank + 1, 1, MPI_COMM_WORLD);
   }


}

