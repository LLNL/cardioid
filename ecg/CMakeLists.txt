
blt_add_executable(
  NAME ecg
  INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
  SOURCES
     ecg.cpp
  DEPENDS_ON mfem simUtil mpi
  )

install(TARGETS ecg
        RUNTIME DESTINATION bin
        )
