#!/bin/bash

if [ "$1" != "" ]; then
	patientDir="${1}"
fi

if [ "${ACVD_DIR}" == "" ]; then
	echo "ACVD_DIR is undefined, please set the ACVD executable directory."
	exit 1
fi

if [ "${TET_DIR}" == "" ]; then
	echo "TET_DIR is undefined, please set the TetGen executable directory."
	exit 1
fi

# grab the trailing patient number from the input directory
patientNum=$(echo "${patientDir}" | grep -o -E "[0-9]+$")

echo "Combining models and clipping mesh..."
python2 vtkClipDataSetWithPolydata.py --path ${patientDir}/
if [ "$?" != 0 ]; then
	echo ""
	echo "Clipping script failed. Aborting..."
	exit 1
fi

echo ""
echo "Fixing mesh..."
echo ""
${ACVD_DIR}/ACVD clipped.vtk 250000 0 -of clipACVD.vtk

if [ "$?" != 0 ]; then
	echo ""
	echo "ACVD failed. Aborting..."
	exit 1
fi

echo ""
echo "Converting VTK to PLY..."
echo ""
${ACVD_DIR}/vtk2ply smooth_clipACVD.vtk

if [ "$?" != 0 ]; then
	echo ""
	echo "VTK to PLY conversion failed. Aborting..."
	exit 1
fi

echo ""
echo "Creating final mesh via TetGen..."
echo ""
${TET_DIR}/tetgen -kCAVNEFq -a0.1666 mesh.ply

if [ "$?" != 0 ]; then
	echo ""
	echo "TetGen failed. Aborting..."
	exit 1
fi

echo ""
echo "Writing output to ${patientDir}/${patientNum}.vtk"
echo ""

if [ -d "${patientDir}" ]; then
	mv mesh.1.vtk ${patientDir}/${patientNum}.vtk
fi
