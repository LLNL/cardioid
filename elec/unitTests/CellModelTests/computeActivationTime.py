#!/usr/bin/python

import sys

def isFloat(str):
    try:
        f = float(str)
        return 1
    except ValueError:
        return 0

def main():
    t1 = 0.0
    v1 = -85.2
    t2 = 0.0
    v2 = -85.2

    for line in sys.stdin:
        cols = line.split()
        if isFloat(cols[0]) == 0:
            continue
        if float(cols[1]) < 0:
            t1 = float(cols[0])
            v1 = float(cols[1])
            continue
        t2 = float(cols[0])
        v2 = float(cols[1])

        activation = t1 +  v1*(t1-t2)/(v1-v2)

        print activation
        break
    




if __name__ == "__main__":
    main()
