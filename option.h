/* 
 * File:   option.h
 * Author: zhang30
 *
 * Created on August 3, 2017, 10:00 AM
 */

#ifndef OPTION_H
#define	OPTION_H

struct Option{
    const char *mesh_file;
    int order;
    bool static_cond;
    bool visualization;
    // run omar's task
    bool omar_task;
    bool omar_fast;
    const char *fiblocs;
    
    // verbose print out
    bool verbose;
    
    // Base angle
    double angle;
    
    double a_endo;
    double a_epi;
    double b_endo;
    double b_epi;
    
    // grid spacing
    double dd;
    // conductivity
    double gL; // mS/mm
    double gT; // mS/mm
    double gN; // mS/mm   

    // cutoff for kdtree point range search rangeCutoff=rcut*maxEdgeLen
    double rcut;  
    double maxEdgeLen;
};

#endif	/* OPTION_H */

