
PROJECT(ViennaSHE)

cmake_minimum_required(VERSION 3.0)

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

   set(PROJECT_SOURCE_DIR_FIXED "${PROJECT_SOURCE_DIR}")  # Looks like CMake behavior has changed over the versions.

   #
   # Full package:
   #
   add_custom_target(dist
      COMMAND "${CMAKE_COMMAND}" -E make_directory "${DISTDIR}"

      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR_FIXED}/LICENSE" "${DISTDIR}/LICENSE"

      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR_FIXED}/changelog" "${DISTDIR}/changelog"

      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR_FIXED}/README" "${DISTDIR}/README"

      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR_FIXED}/CTestConfig.cmake" "${DISTDIR}/CTestConfig.cmake"

      COMMAND "${CMAKE_COMMAND}" -E make_directory "${DISTDIR}/build"


      COMMAND "${CMAKE_COMMAND}" -E make_directory "${DISTDIR}/cmake"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR_FIXED}/cmake/copymanual.cmake" "${DISTDIR}/cmake/copymanual.cmake"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR_FIXED}/cmake/FindOpenCL.cmake" "${DISTDIR}/cmake/FindOpenCL.cmake"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR_FIXED}/cmake/FindPETSc.cmake" "${DISTDIR}/cmake/FindPETSc.cmake"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR_FIXED}/cmake/ViennaSHECommon.cmake" "${DISTDIR}/cmake/ViennaSHECommon.cmake"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR_FIXED}/cmake/ViennaSHEConfig.cmake.in" "${DISTDIR}/cmake/ViennaSHEConfig.cmake.in"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR_FIXED}/cmake/ViennaSHEConfigVersion.cmake.in" "${DISTDIR}/cmake/ViennaSHEConfigVersion.cmake.in"

      COMMAND "${CMAKE_COMMAND}" -E make_directory "${DISTDIR}/doc"
      COMMAND "${CMAKE_COMMAND}" -E copy_directory
         "${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen/html"
         "${DISTDIR}/doc/html"

      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR_FIXED}/doc/Doxyfile.in" "${DISTDIR}/doc/Doxyfile.in"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR_FIXED}/doc/doxyfilter.py" "${DISTDIR}/doc/doxyfilter.py"

      COMMAND "${CMAKE_COMMAND}" -E copy_directory
         "${PROJECT_SOURCE_DIR_FIXED}/src" "${DISTDIR}/src"

      COMMAND "${CMAKE_COMMAND}" -E copy_directory
         "${PROJECT_SOURCE_DIR_FIXED}/libviennashe" "${DISTDIR}/libviennashe"

      COMMAND "${CMAKE_COMMAND}" -E copy_directory
         "${PROJECT_SOURCE_DIR_FIXED}/python" "${DISTDIR}/python"

      COMMAND "${CMAKE_COMMAND}" -E copy_directory
         "${PROJECT_SOURCE_DIR_FIXED}/examples" "${DISTDIR}/examples"

      COMMAND "${CMAKE_COMMAND}" -E make_directory "${DISTDIR}/tests"
      COMMAND "${CMAKE_COMMAND}" -E make_directory "${DISTDIR}/tests/src"
      COMMAND "${CMAKE_COMMAND}" -E copy_directory
         "${PROJECT_SOURCE_DIR_FIXED}/tests/data" "${DISTDIR}/tests/data"
      COMMAND "${CMAKE_COMMAND}" -E copy_directory
         "${PROJECT_SOURCE_DIR_FIXED}/tests/src" "${DISTDIR}/tests/src"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR_FIXED}/tests/CMakeLists.txt" "${DISTDIR}/tests/CMakeLists.txt"

      COMMAND "${CMAKE_COMMAND}" -E copy_directory
         "${PROJECT_SOURCE_DIR_FIXED}/viennacl" "${DISTDIR}/viennacl"

      COMMAND "${CMAKE_COMMAND}" -E copy_directory
         "${PROJECT_SOURCE_DIR_FIXED}/viennagrid" "${DISTDIR}/viennagrid"

      COMMAND "${CMAKE_COMMAND}" -E copy_directory
         "${PROJECT_SOURCE_DIR_FIXED}/viennashe" "${DISTDIR}/viennashe"

      COMMAND "${CMAKE_COMMAND}" -E tar czf ViennaSHE-${VERSION}.tar.gz ViennaSHE-${VERSION}
      WORKING_DIRECTORY "${PROJECT_BINARY_DIR}/dist"
      COMMENT "Generating dist/ViennaSHE-${VERSION}.tar.gz"
      VERBATIM)
   add_dependencies(dist apidoc)

endif(ENABLE_DIST)