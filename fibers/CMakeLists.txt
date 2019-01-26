
blt_add_executable(
  NAME fiberp
  INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
  SOURCES
     io.cpp
     fiberp.cpp
     solver.cpp
     utils.cpp
     triplet.cpp
     genfiber.cpp
     cardfiber.cpp
     cardgradientsp.cpp
  DEPENDS_ON mfem simUtil kdtree mpi
  )

install(TARGETS fiberp
        RUNTIME DESTINATION bin
        )
