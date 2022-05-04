
foreach(device-tutorial half-trigate half-trigate-PETSC mosfet mosfet-PETSC nin-diode-1d nin-diode-2d np-diode-bipolar resistor mosfet-dg nin-diode-2d-hde mos1d-dg )
   add_executable(${device-tutorial} ${device-tutorial}.cpp)
   target_link_libraries(${device-tutorial} shesolvers ${OPENCL_LIBRARIES} ${PETSC_LIBRARIES} ${MPI_mpi_cxx_LIBRARY})

  IF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")  # OpenCL linker settings for MacOS are handled separately
    IF(ENABLE_OPENCL)
      TARGET_LINK_LIBRARIES(${device-tutorial}       OpenCL)
    ENDIF(ENABLE_OPENCL)
  ENDIF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    IF(ENABLE_PETSC)
       target_link_libraries(${device-tutorial}  ${PETSC_LIBRARIES} ${MPI_mpi_cxx_LIBRARY})
    ENDIF(ENABLE_PETSC)
   IF(ENABLE_SERIAL_HDF5)
      TARGET_LINK_LIBRARIES(${device-tutorial} ${HDF5_LIBRARIES})
   ENDIF(ENABLE_SERIAL_HDF5)

endforeach()
