

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)

add_library(ddcMDutil
 	ddcMalloc.c
	error.c
	GridAssignmentObject.c
	hardwareInfo.c
	heap.c
	ioUtils.c
	intQueue.c
	lessThan.c
	match.c
	mpiUtils.c
	object.c
	pio.c
	pioFixedRecordHelper.c
	pioHelper.c
	pioVariableRecordHelper.c
	tagServer.c
	three_algebra.c
	units.c
	utilities.c
)
