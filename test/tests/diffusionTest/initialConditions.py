#!/usr/bin/env python

if __name__=='__main__':
    import re
    import sys

    fff = open(sys.argv[1], "r")

    for ii in range(20):
        line = fff.next()
        m = re.search(r"lrec = (\d+);",line)
        if m:
            lrec = int(m.group(1))
        m = re.search(r"nx = (\d+); ny = (\d+); nz = (\d+);",line)
        if m:
            nx = int(m.group(1))
            ny = int(m.group(2))
            nz = int(m.group(3))
        print line,

    for line in fff:
        m = re.match(r'^\s*(\d+).*?',line)
        if m:
            gid = int(m.group(1))
            ix = gid % nx
            iy = (gid/nx) % ny
            iz = (gid/nx)/ny

            if (ix+iy+iz) % 2:
                out = 1
            else:
                out = 0

            printstring = "%%%dd %%d" % (lrec-3)
            print printstring % (gid,out)
