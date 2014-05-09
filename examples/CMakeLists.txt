
add_subdirectory(devices)
add_subdirectory(libviennashe)

IF(ENABLE_PYTHON_BINDING)
  add_subdirectory(python)
ENDIF(ENABLE_PYTHON_BINDING)

