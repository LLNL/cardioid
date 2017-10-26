#!/usr/bin/env python

import re

def extractNumbersFromLine(line):
    numberReg = re.compile(r'''(?x) 
        [-+]? 
        (?:
           (?:\d+ (?:\.\d*)?) 
           | 
           (?:\.\d+) 
        )
        (?: [eE] [+-]? \d+ )?
    ''')
    numbers = ()
    for match in numberReg.finditer(line):
        numbers += (float(match.group(0)),)
    return numbers

def computeRange(currentMax,currentMin):
    assert(len(currentMax) == len(currentMin))
    currentRange = [None]*len(currentMax)
    for ii in range(0,len(currentMax)):
        currentRange[ii] = currentMax[ii]-currentMin[ii]
    return currentRange

def processErrors(currentDataRun,tolerance):
    globalErrors = []
    currentNumCount = 0

    for ((gl, gln, gn),(bl,bln,bn)) in currentDataRun:
        if not currentNumCount:
            currentNumCount = len(gn)
            currentMax = list(gn)
            currentMin = list(gn)
        for ii in range(0,currentNumCount):
            currentMax[ii] = max(currentMax[ii],gn[ii])
            currentMin[ii] = min(currentMin[ii],gn[ii])

    currentRange = [None]*currentNumCount
    for ii in range(0,currentNumCount):
        currentRange[ii] = currentMax[ii]-currentMin[ii]

    for ((gl, gln, gn),(bl,bln,bn)) in currentDataRun:
        error = [None] * currentNumCount;
        for ii in range(0,currentNumCount):
            denom = currentRange[ii]
            if denom == 0:
                denom = gn[ii]
            if denom == 0:
                if bn[ii] == 0:
                    error[ii] = 0
                else:
                    error[ii] = 1
            else:
                error[ii] = abs(gn[ii]-bn[ii])/denom

        maxError = max(error)
        if maxError > tolerance:
            globalErrors.append("%s | %s | %s\n" %
                                (str(error),
                                 gl.rstrip(),
                                 bl.rstrip()))

    return globalErrors
    
def numCompare(good,bad,tolerance):
    goodLineno = 0
    badLineno = 0
    currentDataRun = []
    currentNumCount = 0
    globalErrors = []
    iterBad = iter(bad)
    doBadEofCheck = False
    for goodLine in good:
        goodNumbers = extractNumbersFromLine(goodLine)
        goodLineno += 1
        if not goodNumbers:
            continue
        
        try:
            foundNumbers = False
            while not foundNumbers:
                badLine = next(iterBad)
                badNumbers = extractNumbersFromLine(badLine)
                badLineno += 1
                foundNumbers = bool(badNumbers)
        except StopIteration:
            globalErrors.append("Premature end-of-file | %s | ---\n" % (goodLine.rstrip()))
            break

        if len(goodNumbers) != len(badNumbers):
            globalErrors.append("Inconsistent columns | %s | %s\n" %
                                (goodLine.rstrip(), badLine.rstrip()))
            break
        if len(goodNumbers) != currentNumCount:
            globalErrors += processErrors(currentDataRun,tolerance)

            #reset the data counters
            currentDataRun = []
            currentNumCount = len(goodNumbers)
        
        currentDataRun.append(((goodLine, goodLineno, goodNumbers),
                               (badLine,badLineno,badNumbers)))

    else:
        doBadEofCheck = True
    if currentDataRun:
        globalErrors += processErrors(currentDataRun,tolerance)
    if doBadEofCheck:
        try:
            badLine = next(iterBad)
            badNumbers = extractNumbersFromLine(badLine)
            badLineno += 1
            if badNumbers:
                globalErrors.append("Premature end-of-file | --- | %s\n" %
                                    badLine)
            
        except StopIteration:
            pass
    return globalErrors


if __name__=='__main__':
    import sys
    import os

    errors = numCompare(open(sys.argv[1],"r"),
                        open(sys.argv[2],"r"),
                        float(sys.argv[3]))
    for line in errors:
        print line,

    sys.exit(len(errors))

    good = ["  Blank line\n",
            " 1  25\n",
            " 2  125\n",
            " 3  27 1e-3 -0.45 +.45",
            " 4  1\n",
            ]
    bad1 = []
    bad2 = [" blank line\n",
            "another blank line\n",
            "1  25.0001\n",
            "2  125.0001\n",
            "following up a blank line"
            " 3  27 1e-3 -0.45 +.45",
            " 4  1.0001\n",
            ]
    bad3 = ["1 "]

    print numCompare(good,bad1, 5e-5)
    print numCompare(good,bad2, 5e-5)
    print numCompare(good,bad2, 5e-7)
    print numCompare(good,bad3, 5e-5)
