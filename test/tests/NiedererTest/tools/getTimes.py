#!/usr/bin/python

import sys

def main():
    
    nx=0
    ny=0
    nz=0
    
    for line in sys.stdin:
#        print line
        if line.find("nx") != -1:
            nx = long(line.split()[2][:-1])
        if line.find("ny") != -1:
            ny = long(line.split()[2][:-1])
        if line.find("nz") != -1:
            nz = long(line.split()[2][:-1])
        if line.find("}") != -1: 
            break


    if nx==0 or ny==0 or nz==0 :
        print "Error finding nx, ny or nz"
        exit(-1)

    xMax = 0
    xMin = 0
    yMax = 0
    yMin = 0
    zMax = 0
    zMin = 0
    P1 = -1.0
    P2 = -1.0
    P3 = -1.0
    P4 = -1.0
    P5 = -1.0
    P6 = -1.0
    P7 = -1.0
    P8 = -1.0
    C  = -1.0
    
    for line in sys.stdin:
        cols = line.split()
        if len(cols) != 2 : continue

        gid = long(cols[0])
        t = float(cols[1])

        x = gid % nx
        gid /= nx
        y = gid % ny
        z = gid / ny

        if x < xMin : xMin = x
        if x > xMax : xMax = x
        if y < yMin : yMin = y
        if y > yMax : yMax = y
        if z < zMin : zMin = z
        if z > zMax : zMax = z

        

        if (x == xMin and y == yMin and z == zMin): P1 = t
        if (x == xMin and y == yMin and z == zMax): P2 = t
        if (x == xMax and y == yMin and z == zMin): P3 = t
        if (x == xMax and y == yMin and z == zMax): P4 = t
        if (x == xMin and y == yMax and z == zMin): P5 = t
        if (x == xMin and y == yMax and z == zMax): P6 = t
        if (x == xMax and y == yMax and z == zMin): P7 = t
        if (x == xMax and y == yMax and z == zMax): P8 = t
        if (x == nx//2    and y == ny//2    and z == nz//2):     C = t

    print "%6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f\\\\"%(P1, P2, P3, P4, P5, P6, P7, P8, C)

if __name__ == "__main__":
    main()
