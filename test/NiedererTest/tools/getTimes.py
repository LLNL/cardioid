#!/usr/bin/python

import sys

def main():
    
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
        if len(cols) != 4 : continue

        x = float(cols[0])
        y = float(cols[1])
        z = float(cols[2])
        t = float(cols[3])

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
        if (x == 0    and y == 0    and z == 0):     C = t

    print "%6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f & %6.2f\\\\"%(P1, P2, P3, P4, P5, P6, P7, P8, C)

if __name__ == "__main__":
    main()
