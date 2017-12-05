#!/usr/bin/env python
import re

class PioData:
    def __init__(self):
        self.header = {}
        self.data = []
        self.foundHeader = False
    def append(self, fileobj):
        inHeader = False
        floatingPoint=r'[+-]?(?:(?:\d+(?:\.\d*)?)|(?:\.\d+))(?:[eE][+-]\d+)?'
        parameterSet = r'(\w+)\s*=\s*((?:\s*?(?:"[^"]")|[^;"]+?)+)\s*;'
        for line in fileobj:
            if inHeader and re.search(parameterSet,line):
                parameters = re.findall(parameterSet, line)
                for (name,value) in parameters:
                    self.header[name] = value
            elif re.search(floatingPoint, line):
                numbers = re.findall(floatingPoint, line)
                gid = int(numbers[0])
                numbers = [float(x) for x in numbers[1:]]
                while len(self.data) < len(numbers):
                    self.data.append({})
                for ii in range(0,len(numbers)):
                    self.data[ii][gid] = numbers[ii]                
            if re.match(r'\w+\s+FILEHEADER\s*\{\s*', line):
                inHeader=True
            if re.match(r'\s*\}\s*', line):
                fieldNames = re.split(r'\s+',self.header["field_names"])
                assert(fieldNames[0] == "gid")
                cursor=0
                self.fieldToId = {}
                for field in fieldNames[1:]:
                    self.fieldToId[field] = cursor
                    cursor += 1
                inHeader=False
                self.foundHeader = True
                
    def getHeader(self, field):
        return self.header[field]
    def getFields(self):
        return self.data.keys()
    def getData(self, gid, field):
        return self.data[self.fieldToId[field]][gid]

def main():
    import sys
    import argparse
    import os
    ap=argparse.ArgumentParser(description="Trace variables from Cardioid output")
    ap.add_argument("--gid", "-g",
                    help="Gid to print",
                    action="append",
                    default=[],
                    )
    ap.add_argument("--field", "-m",
                    help="Field to print",
                    type=str,
                    default="Vm",
                    )
    ap.add_argument("--template", "-t",
                    help="Filename template",
                    type=str,
                    default="Vm",
                    )
    ap.add_argument("directories",
                    nargs=argparse.REMAINDER,
                    help="Directories to search",
                    )
    options = ap.parse_args()
    gids = [int(x) for x in options.gid]
    field = options.field
    template = options.template
    directories = options.directories
    if not directories:
        directories = ["."]
    
    for directoryName in directories:
        #scan this directory for snapshot files
        for snapshotName in [x for x in os.listdir(directoryName)
                             if re.match(r'snapshot\.\d+',x) and os.path.isdir(x)]:
            #scan this directory for files that match the template
            piofile = PioData()
            for datafilename in [x for x in os.listdir(snapshotName)
                                 if re.match(template+r'#\d+',x)]:
                filename=os.path.join(directoryName,snapshotName,datafilename)
                piofile.append(open(filename, "r"))
            if piofile.foundHeader:
                time = float(piofile.getHeader('time'))
                print "%.16g" % time,
                for gid in gids:
                    print "\t%.16g" % piofile.getData(gid,field),
                print ""

if __name__=='__main__':
    main()
