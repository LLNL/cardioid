/* 
 * File:   triplet.cpp
 * Author: zhang30
 *
 * Created on March 2, 2017, 5:28 PM
 */


#include "kdtree++/kdtree.hpp"
#include "triplet.h"

#include <deque>
#include <iostream>
#include <vector>
#include <limits>
#include <functional>
#include <set>

// used to ensure all triplets that are accessed via the operator<< are initialised.
std::set<const void*> registered;

triplet::triplet(value_type a, value_type b, value_type c, int ind) {
    d[0] = a;
    d[1] = b;
    d[2] = c;
    index = ind;
    bool reg_ok = (registered.find(this) == registered.end());
    assert(reg_ok);
    registered.insert(this).second;
}

triplet::triplet(const triplet & x) {
    d[0] = x.d[0];
    d[1] = x.d[1];
    d[2] = x.d[2];
    index = x.index;
    bool reg_ok = (registered.find(this) == registered.end());
    assert(reg_ok);
    registered.insert(this).second;
}

triplet::~triplet() {
    bool unreg_ok = (registered.find(this) != registered.end());
    assert(unreg_ok);
    registered.erase(this);
}

int triplet::getIndex() const{
    return index;
}

double triplet::distance_to(triplet const& x) const {
    double dist = 0;
    for (int i = 0; i != 3; ++i)
        dist += (d[i] - x.d[i])*(d[i] - x.d[i]);
    return std::sqrt(dist);
}

inline bool operator==(triplet const& A, triplet const& B) {
    return A.d[0] == B.d[0] && A.d[1] == B.d[1] && A.d[2] == B.d[2];
}

std::ostream& operator<<(std::ostream& out, triplet const& T) {
    assert(registered.find(&T) != registered.end());
    return out << '(' << T.d[0] << ',' << T.d[1] << ',' << T.d[2] << ", index=" << T.index << ')';
}





