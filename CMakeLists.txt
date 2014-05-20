
cmake_minimum_required(VERSION 2.6)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0009 NEW)
  cmake_policy(SET CMP0011 NEW)
endif(COMMAND cmake_policy)

include(cmake/CMakeLists.txt.in)

add_subdirectory(doc)

# Add visibility of headers
# Necessary for Qt-Creator usage.
################
file( GLOB_RECURSE MAKE_HEADERS_VISIBLE_SRC
    *.hpp *.h *.dox)
add_custom_target( MAKE_HEADERS_VISIBLE SOURCES ${MAKE_HEADERS_VISIBLE_SRC} )

# Dist target
#############

option(ENABLE_DIST "Build the dist package" FALSE)

if(ENABLE_DIST)
   file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/dist")

   set(DISTDIR "${PROJECT_BINARY_DIR}/dist/ViennaSHE-${VERSION}")
   configure_file(cmake/CMakeLists.txt.in "${DISTDIR}/CMakeLists.txt" @ONLY)

   #
   # Full package:
   #
   add_custom_target(dist
      COMMAND "${CMAKE_COMMAND}" -E make_directory "${DISTDIR}"

      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR}/LICENSE" "${DISTDIR}/LICENSE"

      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR}/changelog" "${DISTDIR}/changelog"

      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR}/README" "${DISTDIR}/README"

      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR}/CTestConfig.cmake" "${DISTDIR}/CTestConfig.cmake"

      COMMAND "${CMAKE_COMMAND}" -E make_directory "${DISTDIR}/build"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR}/cmake/README" "${DISTDIR}/build/README"


      COMMAND "${CMAKE_COMMAND}" -E make_directory "${DISTDIR}/cmake"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR}/cmake/copymanual.cmake" "${DISTDIR}/cmake/copymanual.cmake"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR}/cmake/FindOpenCL.cmake" "${DISTDIR}/cmake/FindOpenCL.cmake"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR}/cmake/ViennaSHECommon.cmake" "${DISTDIR}/cmake/ViennaSHECommon.cmake"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR}/cmake/ViennaSHEConfig.cmake.in" "${DISTDIR}/cmake/ViennaSHEConfig.cmake.in"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR}/cmake/ViennaSHEConfigVersion.cmake.in" "${DISTDIR}/cmake/ViennaSHEConfigVersion.cmake.in"

      COMMAND "${CMAKE_COMMAND}" -E make_directory "${DISTDIR}/doc"
      COMMAND "${CMAKE_COMMAND}" -E copy_directory
         "${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen/html"
         "${DISTDIR}/doc/html"

      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR}/doc/Doxyfile.in" "${DISTDIR}/doc/Doxyfile.in"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR}/doc/doxyfilter.py" "${DISTDIR}/doc/doxyfilter.py"

      COMMAND "${CMAKE_COMMAND}" -E copy_directory
         "${PROJECT_SOURCE_DIR}/src" "${DISTDIR}/src"

      COMMAND "${CMAKE_COMMAND}" -E copy_directory
         "${PROJECT_SOURCE_DIR}/libviennashe" "${DISTDIR}/libviennashe"

      COMMAND "${CMAKE_COMMAND}" -E copy_directory
         "${PROJECT_SOURCE_DIR}/python" "${DISTDIR}/python"

      COMMAND "${CMAKE_COMMAND}" -E copy_directory
         "${PROJECT_SOURCE_DIR}/examples" "${DISTDIR}/examples"

      COMMAND "${CMAKE_COMMAND}" -E copy_directory
         "${PROJECT_SOURCE_DIR}/external" "${DISTDIR}/external"

      COMMAND "${CMAKE_COMMAND}" -E make_directory "${DISTDIR}/tests"
      COMMAND "${CMAKE_COMMAND}" -E make_directory "${DISTDIR}/tests/src"
      COMMAND "${CMAKE_COMMAND}" -E copy_directory
         "${PROJECT_SOURCE_DIR}/tests/data" "${DISTDIR}/tests/data"
      COMMAND "${CMAKE_COMMAND}" -E copy_directory
         "${PROJECT_SOURCE_DIR}/tests/src" "${DISTDIR}/tests/src"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR}/tests/CMakeLists.txt" "${DISTDIR}/tests/CMakeLists.txt"

      COMMAND "${CMAKE_COMMAND}" -E copy_directory
         "${PROJECT_SOURCE_DIR}/viennacl" "${DISTDIR}/viennacl"

      COMMAND "${CMAKE_COMMAND}" -E copy_directory
         "${PROJECT_SOURCE_DIR}/viennagrid" "${DISTDIR}/viennagrid"

      COMMAND "${CMAKE_COMMAND}" -E copy_directory
         "${PROJECT_SOURCE_DIR}/viennashe" "${DISTDIR}/viennashe"

      COMMAND "${CMAKE_COMMAND}" -E tar czf ViennaSHE-${VERSION}.tar.gz ViennaSHE-${VERSION}
      WORKING_DIRECTORY "${PROJECT_BINARY_DIR}/dist"
      COMMENT "Generating dist/ViennaSHE-${VERSION}.tar.gz"
      VERBATIM)
   add_dependencies(dist apidoc manual)

endif(ENABLE_DIST)
