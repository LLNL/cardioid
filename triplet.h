/* 
 * File:   triplet.h
 * Author: zhang30
 *
 * Created on March 2, 2017, 5:28 PM
 */

#ifndef TRIPLET_H
#define	TRIPLET_H

#include "kdtree++/kdtree.hpp"

#include <deque>
#include <iostream>
#include <vector>
#include <limits>
#include <functional>
#include <set>

struct triplet {
    typedef double value_type;

    triplet(value_type a, value_type b, value_type c, int ind);

    triplet(const triplet & x);

    ~triplet();

    double distance_to(triplet const& x) const;
    
    int getIndex() const;

    inline value_type operator[](size_t const N) const {
        return d[N];
    };

    value_type d[3];
    int index;
};

inline bool operator==(triplet const& A, triplet const& B);
std::ostream& operator<<(std::ostream& out, triplet const& T);

inline double tac(triplet t, size_t k) {
    return t[k];
}

typedef KDTree::KDTree<3, triplet, std::pointer_to_binary_function<triplet,size_t,double> > tree_type;

#endif	/* TRIPLET_H */

