
blt_add_executable(
  NAME mechanics_driver
  INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
  SOURCES
     ConstantTension.cpp
     Lumens2009.cpp
     ExcitationContraction.cpp
     cardiac_coefficients.cpp
     cardiac_integrators.cpp
     cardiac_physics.cpp
     cardiac_solvers.cpp
     mechanics_driver.cpp
  DEPENDS_ON mfem mpi
  )
