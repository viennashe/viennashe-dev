
add_library(viennashe SHARED src/libviennashe.cpp src/device.cpp src/material.cpp src/config.cpp src/simulator.cpp src/output.cpp src/quantity.cpp )

IF(MSVC)
set_source_files_properties(src/simulator.cpp PROPERTIES COMPILE_FLAGS /bigobj)
ENDIF(MSVC)


target_link_libraries(viennashe shesolvers)
target_link_libraries(viennashe ${OPENCL_LIBRARIES} ${PETSC_LIBRARIES} ${MPI_mpi_cxx_LIBRARY})

