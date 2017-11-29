#!/usr/bin/env python

import random
import math
import numpy

def quat2rot(q):
    w = q[0]
    x = q[1]
    y = q[2]
    z = q[3]

    x2=x*x;
    y2=y*y;
    z2=z*z;
    xy=x*y;
    xz=x*z;
    yz=y*z;
    wx=w*x;
    wy=w*y;
    wz=w*z;

    V = numpy.zeros((3,3))
    V[0,0] = 1-2*y2-2*z2
    V[1,0]=2*xy-2*wz;
    V[2,0]=2*xz+2*wy;
    
    V[0,1]=2*xy+2*wz;
    V[1,1]=1-2*x2-2*z2;
    V[2,1]=2*yz-2*wx;
    
    V[0,2]=2*xz-2*wy;
    V[1,2]=2*yz+2*wx;
    V[2,2]=1-2*x2-2*y2;

    return V


def main():
    import sys
    import argparse
    ap=argparse.ArgumentParser(description="Make a mesh with missing voxels to test diffusion")
    ap.add_argument("--size", "-n",
                    help="Size of the problem",
                    type=int,
                    default=40,
                    )
    ap.add_argument("--gil",
                    help="longitudinal conductivity",
                    type=float,
                    default=1.0,
                    )
    ap.add_argument("--git",
                    help="transverse conductivity",
                    type=float,
                    default=0.1,
                    )
    ap.add_argument("--gin",
                    help="normal conductivity",
                    type=float,
                    default=0.01,
                    )
    ap.add_argument("--seed", "-z",
                    help="Initial random seed",
                    type=int,
                    default=None,
                    )
    
    options = ap.parse_args()
    size=options.size
    random.seed(options.seed)
    gil = options.gil
    git = options.git
    gin = options.gin
    #outfilename = options.outfilename

    usedGids = set()
    for iz in range(size):
        for iy in range(size):
            for ix in range(size):
                gid = ix+size*(iy+size*iz)
                probability = float(ix+iy+iz)/(size-1)/3
                if random.random() > probability:
                    usedGids.add(gid)
    
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
""" % (len(usedGids),size,size,size)

    for gid in usedGids:
        quaternion = [0,0,0,0]
        maxLength=math.sqrt(3)
        quaternion[0] = (2*random.random()-1)/maxLength
        quaternion[1] = (2*random.random()-1)/maxLength
        quaternion[2] = (2*random.random()-1)/maxLength
        quaternion[3] = math.sqrt(1-quaternion[0]**2-quaternion[1]**2-quaternion[2]**2)
        V = quat2rot(quaternion)
        cond = V.dot(numpy.diag([gil, git, gin]).dot(V.transpose()))
        print "%10d 101 %21.16g %21.16g %21.16g %21.16g %21.16g %21.16g" % (
            gid, cond[0,0], cond[0,1], cond[0,2], cond[1,1], cond[1,2], cond[2,2],
            )

if __name__=='__main__':
    main()
