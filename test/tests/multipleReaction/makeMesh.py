#!/usr/bin/env python

import random
import math
import numpy
import bisect

def main():
    import sys
    import argparse
    ap=argparse.ArgumentParser(description="Make a slab with 3 tissue types")
    ap.add_argument("--nx", "-x",
                    help="Size of the problem",
                    type=int,
                    default=40,
                    )
    ap.add_argument("--ny", "-y",
                    help="Size of the problem",
                    type=int,
                    default=40,
                    )
    ap.add_argument("--nz", "-z",
                    help="Size of the problem",
                    type=int,
                    default=40,
                    )
    ap.add_argument("--gil",
                    help="longitudinal conductivity",
                    type=float,
                    default=0.1334177,
                    )
    ap.add_argument("--git",
                    help="transverse conductivity",
                    type=float,
                    default=0.0176062,
                    )
    ap.add_argument("--gin",
                    help="normal conductivity",
                    type=float,
                    default=0.0176062,
                    )
    
    options = ap.parse_args()
    nx=options.nx
    ny=options.ny
    nz=options.nz
    gil = options.gil
    git = options.git
    gin = options.gin
    tissueTypes = [100,101,102]
    #outfilename = options.outfilename
    
    print """anatomy FILEHEADER {
  datatype = VARRECORDASCII;
  nfiles = 1;
  nrecords = %d;
  nfields = 8;
  field_names = gid cellType sigma11 sigma12 sigma13 sigma22 sigma23 sigma33;
  field_types = u u f f f f f f;
  nx = %d; ny = %d; nz = %d;
  field_units = 1 1 mS/mm mS/mm mS/mm mS/mm mS/mm mS/mm;
}
""" % (nx*ny*nz,nx,ny,nz)

    extents = [0]
    for itype in range(0,len(tissueTypes)):
        mySize = (nx/len(tissueTypes))
        if itype < (nx % len(tissueTypes)):
            mySize += 1
        extents.append(extents[-1]+mySize)
    assert(extents[-1] == nx)
        
    gid = 0
    for iz in range(0,nz):
        for iy in range(0,ny):
            for ix in range(0,nx):
                index = bisect.bisect_right(extents,ix)-1
                print "%10d %3d %21.16g 0 0 %21.16g 0 %21.16g" % (
                    gid, tissueTypes[index], gil, git, gin
                )
                
                gid += 1

    sensor = open("sensor.txt", "w")
    for itype in range(0,len(tissueTypes)):
        ix = (extents[itype]+extents[itype+1])/2
        iy = ny/2
        iz = nz/2
        gid = ix+nx*(iy+ny*iz)
        print >>sensor, gid
        

                
if __name__=='__main__':
    main()
