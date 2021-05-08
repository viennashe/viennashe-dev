
add_library(shesolvers solver.cpp)
#message(${PETSC_INCLUDES})
target_include_directories(shesolvers SYSTEM PUBLIC ${PETSC_INCLUDES})
target_include_directories(shesolvers SYSTEM PUBLIC ${MPI_CXX_HEADER_DIR} ${MPI_C_HEADER_DIR})

#target_link_libraries(shesolvers petsc)
IF(NOT MSVC)
  set_source_files_properties(solver.cpp PROPERTIES COMPILE_FLAGS -fPIC)
ENDIF(NOT MSVC)

