import vtk
import numpy as np
import sys
import argparse
import os

# Create polydata to slice the grid with. In this case, use a cone. This could
# be any polydata including a stl file.
# cone = vtk.vtkConeSource()
# cone.SetResolution(20)
# cone.SetCenter(50, 50, 50)
# cone.SetRadius(25)
# cone.SetHeight(40)
# cone.Update()

def distanceForModel(filename):
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(filename)
    reader.Update()

    # implicit function that will be used to slice the mesh
    distance = vtk.vtkImplicitPolyDataDistance()
    distance.SetInput(reader.GetOutput())

    return distance

def clipData(inputData,distanceFunc,insideOut):
    evalDistanceOnGrid("SignedDistance",inputData,distanceFunc)
    clipper = vtk.vtkClipDataSet()
    clipper.SetInputData(inputData)
    if insideOut:
        clipper.InsideOutOn()
    clipper.SetValue(0.0)
    clipper.Update()
    return clipper.GetOutput()

def evalDistanceOnGrid(name,rgrid,distance):
    # Create an array to hold distance information
    signedDistances = vtk.vtkFloatArray()
    signedDistances.SetNumberOfComponents(1)
    signedDistances.SetName(name)

    # Evaluate the signed distance function at all of the grid points
    for pointId in range(rgrid.GetNumberOfPoints()):
        p = rgrid.GetPoint(pointId)
        signedDistance = distance.EvaluateFunction(p)
        signedDistances.InsertNextValue(signedDistance)
        # print(p, signedDistance)

        # add the SignedDistances to the grid
        rgrid.GetPointData().SetScalars(signedDistances)

# Note: Since the RV and LV are contained within the epicardium, we only need to find min/max of the epicardium for use in the bounding box
# This function just reads the epi fiducial file and returns a list of the mins and maxes for each direction
def bBoxBounds(epi):
    try:
        with open(epi) as file:
            mins = [1000] * 3
            maxs = [-1000] * 3
            for line in file:
                lineData = line.split()
                if lineData[0] != "#":
                    lineData = line.split(",")
                    for i in range(1, 4): # 3D
                        temp = float(lineData[i])
                        if temp < mins[i - 1]:
                            mins[i - 1] = temp
                        if temp > maxs[i - 1]:
                            maxs[i - 1] = temp

        return mins, maxs

    except IOError as err:
        #print err
        sys.exit()

# Computes the direction and value of the desired cut
def cuttingPlane(cut):
    # can this be done more efficiently?
    success = False
    try:
        with open(cut) as file:
            #print(f"Opening {cut} ...")
            success = True
            pts = []
            for line in file:
                lineData = line.split()
                if lineData[0] != "#":
                    lineData = line.split(",")
                    pts.append(lineData[1:4])

    except Exception as err:
        #print err
        sys.exit()

    if success:
        diff = []
        for i in range(len(pts)):
            diffTemp = [0x0, 0x0, 0x0]
            for j in range(len(pts[0])):
                pts[i][j] = float(pts[i][j])
                # don't compute difference for the first row
                if i > 0:
                    diffTemp[j] = pts[i][j] - pts[i - 1][j]

            if i > 0:
                diff.append(diffTemp)

        # Go through the diff list and compute the direction that has the smallest average difference, and use that as the plane direction. Also, use the average value as the plane bound.
        #avg = np.mean(diff)
        avg = []
        for i in range(len(diff[0])):
            avg.append(np.mean([abs(col[i]) for col in diff]))

        # find the minimum of the average differences, and extract its corresponding index
        val, planeDirection = min((val, idx) for (idx, val) in enumerate(avg))

        # use the mean of the coordinate specified by planeDirection as the planeBound
        planeBound = np.mean([col[planeDirection] for col in pts])

        return planeDirection, planeBound




parser = argparse.ArgumentParser()

parser.add_argument("-p", "--path", help = "path for all the necessary data", type = str)
args = parser.parse_args()
pathname = args.path

fiducialLV = os.path.join(pathname, "LVF.fcsv")
fiducialRV = os.path.join(pathname, "RVF.fcsv")
fiducialEPI = os.path.join(pathname, "EF.fcsv")
fiducialCut = os.path.join(pathname, "cutting.fcsv")

vtkLV = os.path.join(pathname, "lv.vtk")
vtkRV = os.path.join(pathname, "rv.vtk")
vtkEPI = os.path.join(pathname, "epi.vtk")

# fiducialLV = "lv/LVF.fcsv"
# fiducialRV = "rv/RVF.fcsv"
# fiducialEPI = "epi/EF.fcsv"
# fiducialCut = "cutting.fcsv"

lvDistance = distanceForModel(vtkLV)
rvDistance = distanceForModel(vtkRV)
epiDistance = distanceForModel(vtkEPI)

mins, maxs = bBoxBounds(fiducialEPI)
direction, bound = cuttingPlane(fiducialCut)

diffMin = abs(mins[direction] - bound)
diffMax = abs(maxs[direction] - bound)

if diffMin < diffMax:
    bboxIndex = 0
else:
    bboxIndex = 1

buffer = 0.05 # 5% buffer
# create a grid
boundingBox = np.array([[mins[0], maxs[0]],
                        [mins[1], maxs[1]],
                        [mins[2], maxs[2]]])

for i in range(len(boundingBox)):
    boundingBox[i] += boundingBox[i] * buffer

boundingBox[direction][bboxIndex] = bound

# print mins, maxs

# boundingBox = np.array([[-140,-55],
#                         [ -30, 80],
#                         [-120, 10]]);

boxSize = boundingBox[:,1]-boundingBox[:,0]
volume = float(np.prod(boxSize))

numPoints = 2000000

h = (volume/numPoints)**(1./3)

xCoords = vtk.vtkFloatArray()
for x, i in enumerate(np.linspace(boundingBox[0,0], boundingBox[0,1], round(boxSize[0]//h))):
    xCoords.InsertNextValue(i)

yCoords = vtk.vtkFloatArray()
for y, i in enumerate(np.linspace(boundingBox[1,0], boundingBox[1,1], round(boxSize[1]//h))):
    yCoords.InsertNextValue(i)

zCoords = vtk.vtkFloatArray()
for z, i in enumerate(np.linspace(boundingBox[2,0], boundingBox[2,1], round(boxSize[2]//h))):
    zCoords.InsertNextValue(i)

# The coordinates are assigned to the rectilinear grid. Make sure that
# the number of values in each of the XCoordinates, YCoordinates,
# and ZCoordinates is equal to what is defined in SetDimensions().

rgrid = vtk.vtkRectilinearGrid()
rgrid.SetDimensions(x+1, y+1, z+1)
rgrid.SetXCoordinates(xCoords)
rgrid.SetYCoordinates(yCoords)
rgrid.SetZCoordinates(zCoords)



#evalDistanceOnGrid("lv",rgrid,lvDistance)
#evalDistanceOnGrid("rv",rgrid,lvDistance)
#evalDistanceOnGrid("epi",rgrid,lvDistance)

# use vtkClipDataSet to slice the grid with the polydata



clippedOutput = clipData(rgrid,epiDistance,True)
clippedOutput = clipData(clippedOutput,lvDistance,False)
clippedOutput = clipData(clippedOutput,rvDistance,False)

surfaceFilter = vtk.vtkDataSetSurfaceFilter()
surfaceFilter.SetInputData(clippedOutput)
surfaceFilter.Update()
clippedOutput = surfaceFilter.GetOutput()

#print clippedOutput

# --- mappers, actors, render, etc. ---
# mapper and actor to view the cone
#coneMapper = vtk.vtkPolyDataMapper()
#coneMapper.SetInputConnection(reader.GetOutputPort())
#coneActor = vtk.vtkActor()
#coneActor.SetMapper(coneMapper)

# geometry filter to view the background grid
geometryFilter = vtk.vtkRectilinearGridGeometryFilter()
geometryFilter.SetInputData(rgrid)
geometryFilter.SetExtent(0, x+1, 0, y+1, (z+1)//2, (z+1)//2)
geometryFilter.Update()

rgridMapper = vtk.vtkPolyDataMapper()
rgridMapper.SetInputConnection(geometryFilter.GetOutputPort())

wireActor = vtk.vtkActor()
wireActor.SetMapper(rgridMapper)
wireActor.GetProperty().SetRepresentationToWireframe()
wireActor.GetProperty().SetColor(.1, .1, .1)

# mapper and actor to view the clipped mesh
clipperMapper = vtk.vtkDataSetMapper()
clipperMapper.SetInputData(clippedOutput)

clipperActor = vtk.vtkActor()
clipperActor.SetMapper(clipperMapper)
clipperActor.GetProperty().SetRepresentationToWireframe()
clipperActor.GetProperty().SetColor(.1, .1, .1)

# A renderer and render window
renderer = vtk.vtkRenderer()
renderer.SetBackground(1, 1, 1)

# add the actors
#renderer.AddActor(coneActor)
renderer.AddActor(wireActor)
renderer.AddActor(clipperActor)

renwin = vtk.vtkRenderWindow()
renwin.AddRenderer(renderer)

# An interactor
interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(renwin)

# Start
# interactor.Initialize()
# interactor.Start()

writer = vtk.vtkPolyDataWriter()
writer.SetInputData(clippedOutput)
writer.SetFileName("clipped.vtk")
writer.Update()
