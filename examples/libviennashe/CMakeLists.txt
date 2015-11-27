
foreach(device-tutorial resistor1d-c  )
   add_executable(${device-tutorial} ${device-tutorial}.c)
   target_link_libraries(${device-tutorial} viennashe viennagrid)

endforeach()

