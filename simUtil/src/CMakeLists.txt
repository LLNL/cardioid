
blt_add_library(
    NAME simUtil
    DEFINES -DDiff_Weight_Type_Single -DWITH_PIO -DWITH_MPI
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/../include
    SOURCES 
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
    DEPENDS_ON mpi
)
