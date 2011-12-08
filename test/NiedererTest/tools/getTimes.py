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

    for line in sys.stdin:
        cols = line.split()
        if len(cols) != 4 : continue
        if float(cols[0]) < xMin : xMin = float(cols[0])
        if float(cols[0]) > xMax : xMax = float(cols[0])
        if float(cols[1]) < yMin : yMin = float(cols[1])
        if float(cols[1]) > yMax : yMax = float(cols[1])
        if float(cols[2]) < zMin : zMin = float(cols[2])
        if float(cols[2]) > zMax : zMax = float(cols[2])

        if (float(cols[0]) == xMin and float(cols[1]) == yMin and float(cols[2]) == zMin): P1 = float(cols[3])
        if (float(cols[0]) == xMin and float(cols[1]) == yMin and float(cols[2]) == zMax): P2 = float(cols[3])
        if (float(cols[0]) == xMax and float(cols[1]) == yMin and float(cols[2]) == zMin): P3 = float(cols[3])
        if (float(cols[0]) == xMax and float(cols[1]) == yMin and float(cols[2]) == zMax): P4 = float(cols[3])
        if (float(cols[0]) == xMin and float(cols[1]) == yMax and float(cols[2]) == zMin): P5 = float(cols[3])
        if (float(cols[0]) == xMin and float(cols[1]) == yMax and float(cols[2]) == zMax): P6 = float(cols[3])
        if (float(cols[0]) == xMax and float(cols[1]) == yMax and float(cols[2]) == zMin): P7 = float(cols[3])
        if (float(cols[0]) == xMax and float(cols[1]) == yMax and float(cols[2]) == zMax): P8 = float(cols[3])

    print (P1, P2, P3, P4, P5, P6, P7, P8)

if __name__ == "__main__":
    main()
