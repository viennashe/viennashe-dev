
file(GLOB PythonExampleFiles "${CMAKE_CURRENT_SOURCE_DIR}/*.py" )

add_custom_target(python-examples ALL)
foreach(PyExampleFile ${PythonExampleFiles})
  GET_FILENAME_COMPONENT(PyFilename "${PyExampleFile}" NAME )
  ##message("Copying ${PyExampleFile} TO ${CMAKE_CURRENT_BINARY_DIR}/${PyFilename} " )
  add_custom_command(TARGET python-examples PRE_BUILD COMMAND ${CMAKE_COMMAND} -E copy "${PyExampleFile}" "${CMAKE_CURRENT_BINARY_DIR}")
endforeach()

