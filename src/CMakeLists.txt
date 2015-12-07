
# Solver module
add_library(shesolvers solver.cpp)
IF(NOT MSVC)
  set_source_files_properties(solver.cpp PROPERTIES COMPILE_FLAGS -fPIC)
ENDIF(NOT MSVC)


# Boxmesh module
add_library(sheboxmesh boxmesh.cpp)
IF(NOT MSVC)
  set_source_files_properties(boxmesh.cpp PROPERTIES COMPILE_FLAGS -fPIC)
ENDIF(NOT MSVC)
