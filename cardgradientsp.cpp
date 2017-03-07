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
        triplet vetexNearPt = *found.first;
        int vertex = vetexNearPt.getIndex();
        vector<int> elements = vert2Elements[vertex];
        for (unsigned e = 0; e < elements.size(); e++) {
            int eleIndex = elements[e];
            //For barycentric
            Vector q(4);
            q(0) = x;
            q(1) = y;
            q(2) = z;
            q(3) = 1.0;
            if (isInTetElement(q, mesh, eleIndex)) {

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

                DenseMatrix Sigma(dim3, dim3);
                calcSigma(Sigma, QPfib, conduct);

                double *s = Sigma.Data();
                anatomy anat;
                anat.gid = g;
                anat.celltype = getCellType(phi_epi, phi_lv, phi_rv);
                anat.sigma[0] = s[0];
                anat.sigma[1] = s[3];
                anat.sigma[2] = s[6];
                anat.sigma[3] = s[4];
                anat.sigma[4] = s[7];
                anat.sigma[5] = s[8];
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
    int nFiles_ = num_procs;
//    if (nFiles_ > 0)
//        PioSet(file, "ngroup", nFiles_);  

    if (myid == 0) {
        Pprintf(file, "anatomy FILEHEADER { \n");
        Pprintf(file, "  datatype = VARRECORDASCII;\n");
        Pprintf(file, "  nfiles = %d;  \n", nFiles_);
        Pprintf(file, "  nrecord = %d; \n", header.nrecord);
        Pprintf(file, "  nfields = 1; \n");
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

