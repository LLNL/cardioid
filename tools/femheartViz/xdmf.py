
import numpy as np
import h5py
import re
import sys

class Indenter:
    def __init__(self, fff=sys.stdout, indentString="   "):
        self.indent = indentString
        self.indentAmount = 0
        self.outfile = fff

    def __call__(self, string, *args, **kwargs):
        initialstrip = string.lstrip('\n')
        if not initialstrip:
            print(string, file=self.outfile)
            return
        indentString = self.indent*self.indentAmount
        for line in initialstrip.split('\n'):
            outline = indentString+line
            if line == "":
                print("\n", file=self.outfile)
            elif kwargs:
                print((outline % kwargs), file=self.outfile)
            elif len(args) == 1 and type(args[0]) == dict:
                print((outline % args[0]), file=self.outfile)
            elif len(args) == 1:
                print((outline % args[0]), file=self.outfile)
            elif len(args) == 0:
                print(outline, file=self.outfile)
            else:
                print((outline % args), file=self.outfile)

    def inc(self, indentAmount=1):
        self.indentAmount += indentAmount

    def dec(self, indentAmount=1):
        self.indentAmount -= indentAmount

def order(inSearchable):
    aa = [ii for ii in inSearchable]
    aa.sort()
    for ii in aa:
        yield ii
        
        
def outputGrid(out, geom, data, time=None):
    out('<Grid GridType="Uniform">')
    out.inc()
    if time!=None:
        out('<Time Value="{}" />'.format(time))
    out(geom)
    out(data)
    out.dec()
    out('</Grid>')

def attrDesc(h5array, numPoints, numCells):
    attr = ""
    if h5array.shape[0] == numPoints:
        attr += 'Center="Node"'
    elif h5array.shape[0] == numCells:
        attr += 'Center="Cell"'
    else:
        return ""
    if len(h5array.shape) == 1 or h5array.shape[-1] == 1:
        attr += ' AttributeType="Scalar"'
    elif h5array.shape[-1] == 3:
        attr += ' AttributeType="Vector"'
    elif h5array.shape[-1] == 9:
        attr += ' AttributeType="Tensor"'
    elif h5array.shape[-1] == 6:
        attr += ' AttributeType="Tensor6"'
    else:
        return ""
    return attr

def actualTimeFromDir(name):
    strTime = re.sub('^tm0*','', name)
    if not strTime:
        strTime="0"
    return float(strTime)/1000.

def hdf5ArrayDesc(h5array):
    if h5array.dtype.kind == 'i':
        ttype = "Int"
    else:
        ttype = "Float"
    return """   <DataItem Format="HDF" Dimensions="{dim}" NumberType="{ttype}" Precision="{prec}">
      {filename}:{name}
   </DataItem>""".format(dim=" ".join([str(x) for x in h5array.shape]),
           ttype=ttype,
           prec=h5array.dtype.itemsize,
           filename=h5array.file.filename,
           name=h5array.name,
           )

def isH5Array(h5array):
    return re.search("Dataset", str(h5array.__class__))

def main():
    import argparse
    import os
    parser = argparse.ArgumentParser(description="""

This script makes an hdf5 file suitable for reading in paraview or
visit onto stdout.  To work the command needs:

- At most one HDF5 file that contains an array named "XYZ" that lists
  the points in the mesh. The array should be number_of_points x 3.

- At most one HDF5 file that contains an array called "tets" or
  "elements".  If the array is "tets", the array should be integer of
  size number_of_cells x 4.  If the array is "elements", the array
  must be one dimensional with a mixed list of the elements in the
  XDMF format. The elements array must also have an attribute called
  "numCells" that contains the number of cells in the array.

- Any HDF5 files with (1,3,6,9) x (numPoints, numCells) dimensional
  arrays. These arrays will appear as point or cell centered vertices
  respectively.

- Any HDF5 files with "tm\d+" or "tm\d+.\d+", followed by arrays.
  These will show up as time arrays and will make an hdf5 file that
  can be animated.
""")
    parser.add_argument("h5files", nargs='+',
                        help="h5 files to combine into an xdmf")
    args = parser.parse_args()

    geom = ""
    numPoints=-1
    numCells=-1
    for h5filename in args.h5files:
        h5file = h5py.File(h5filename, 'r')
        if 'XYZ' in h5file:
            numPoints = h5file['XYZ'].shape[0]
            geom += """<Geometry GeometryType="XYZ">
{hdf}
</Geometry>
""".format(hdf=hdf5ArrayDesc(h5file['XYZ']))
        if 'tets' in h5file:
            numCells = h5file['tets'].shape[0]
            geom += """<Topology TopologyType="Tetrahedron" NumberOfElements="{numCells}">
{hdf}
</Topology>
""".format(numCells=numCells, hdf=hdf5ArrayDesc(h5file['tets']))
        elif 'elements' in h5file:
            numCells = h5file['elements'].attrs["numCells"]
            geom += """<Topology TopologyType="Mixed" NumberOfElements="{numCells}">
{hdf}
</Topology>
""".format(numCells=numCells, hdf=hdf5ArrayDesc(h5file['elements']))
        h5file.close()

    constData = ""
    timeData = {}
    for h5filename in args.h5files:
        h5file = h5py.File(h5filename, 'r')
        for key in h5file.keys():
            if key in ('XYZ', 'tets', 'elements'):
                continue
            if isH5Array(h5file[key]):
                attr = attrDesc(h5file[key],numPoints,numCells)
                if attr:
                    constData += """<Attribute Name="{key}" {attr}>
{hdf}
</Attribute>
""".format(key=key,
           hdf=hdf5ArrayDesc(h5file[key]),
           attr=attr,)
            elif re.match("^tm[0-9]+(?:\.[0-9]+)?$",key):
                ttt = key
                for keyAtTime in h5file[ttt]:
                    if isH5Array(h5file[ttt][keyAtTime]):
                        attr = attrDesc(h5file[ttt][keyAtTime],numPoints,numCells)
                        if attr:
                            if ttt not in timeData:
                                timeData[ttt] = ""
                            timeData[ttt] += """<Attribute Name="{keyAtTime}" {attr}>
{hdf}
</Attribute>
""".format(keyAtTime=keyAtTime,
           hdf=hdf5ArrayDesc(h5file[ttt][keyAtTime]),
           attr=attr,)
    
    out = Indenter()
    out('''<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">
   <Domain>''')
    out.inc()
    out.inc()
    if not timeData:
        outputGrid(out, geom, constData)
    else:
        out('<Grid GridType="Collection" CollectionType="Temporal">')
        out.inc()
        for (time,timeData) in order(timeData.items()):
            outputGrid(out, geom, constData+timeData, time=actualTimeFromDir(time))
        out.dec()
        out('</Grid>')
    out.dec()
    out("</Domain>")
    out.dec()
    out("</Xdmf>")

if __name__=='__main__':
    main()
