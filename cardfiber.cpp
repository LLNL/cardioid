#define KDTREE_DEFINE_OSTREAM_OPERATORS   // Position sensitive.

#include "cardfiber.h"
#include "constants.h"
#include "triplet.h"
#include "kdtree++/kdtree.hpp"

void buildKDTree(Mesh *mesh, tree_type& kdtree, vector<Vector>& boundingbox, double dd) {
    int NumOfVertices = mesh->GetNV();

    //for (int i = 0; i < NumOfVertices; i++)
    for (int i = 0; i < 10; i++) {
        const double* coord = mesh->GetVertex(i);
        triplet node(coord[0], coord[1], coord[2], i);
        kdtree.insert(node);
    }
    cout << kdtree << endl;
    for (tree_type::const_iterator target = kdtree.begin(); target != kdtree.end(); ++target) {
        cout << *target << endl;
    }
}


