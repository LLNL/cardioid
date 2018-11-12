#ifndef WRITE_CELLS_HH
#define WRITE_CELLS_HH

#include <string>
#include <vector>
class Simulate;
class AnatomyCell;

// This signature is for use during load balancing, when the Simulate
// and Anatomy objects are both poorly formed.  In particular, the
// nLocal() method of Anatomy can't be relied on to give the correct
// number.  Bad things will happen if you call this signature later.
// Once the Simulation is fully set up with the halo exchange the cell
// array has both local and remote cells and this routine will write all
// of them into the file producing a file with duplicate cells.
void writeCells(const std::vector<AnatomyCell>& cells,
                int nx, int ny, int nz,
                const std::string& filename);

#endif
