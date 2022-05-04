

include(CTest)
include(CMakeDependentOption)

# Installation directories
##########################

set(INSTALL_INCLUDE_DIR include CACHE PATH
   "Installation directory for headers")
if(WIN32 AND NOT CYGWIN)
   set(DEF_INSTALL_CMAKE_DIR CMake)
else()
   set(DEF_INSTALL_CMAKE_DIR lib/cmake/viennashe)
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH
   "Installation directory for CMake files")

if(NOT IS_ABSOLUTE "${INSTALL_CMAKE_DIR}")
   set(INSTALL_CMAKE_DIR "${CMAKE_INSTALL_PREFIX}/${INSTALL_CMAKE_DIR}")
endif()
file(RELATIVE_PATH CONF_REL_INSTALL_PREFIX "${INSTALL_CMAKE_DIR}"
   "${CMAKE_INSTALL_PREFIX}")
if(NOT IS_ABSOLUTE "${INSTALL_INCLUDE_DIR}")
   set(INSTALL_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/${INSTALL_INCLUDE_DIR}")
endif()
file(RELATIVE_PATH CONF_REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}"
   "${INSTALL_INCLUDE_DIR}")

# User options
##############

OPTION(BUILD_TESTING "Build the tests " ON)

##OPTION(WITH_VSHE "Build the main ViennaSHE application" ON)

OPTION(ENABLE_PYTHON_BINDINGS "Enable Python bindings. Requires SWIG 2.0" OFF)

OPTION(ENABLE_OPENCL "Enable OpenCL-accelerated solver" OFF)

OPTION(ENABLE_OPENMP "Enable OpenMP-accelerated solver" OFF)

# If you want to build the examples that use Eigen
option(DISABLE_LOGGING "Disables all logging in ViennaSHE" OFF)

## USE THIS TO INCLUDE HDF5 !
##option(ENABLE_SERIAL_HDF5 "Enables serial (no mpi) HDF5 in ViennaSHE" OFF)

mark_as_advanced(DISABLE_LOGGING)


# Find prerequisites
####################

if (ENABLE_OPENCL)
  INCLUDE_DIRECTORIES("external/")
  find_package(OpenCL)
  ADD_DEFINITIONS( -DVIENNACL_WITH_OPENCL )
endif(ENABLE_OPENCL)

if (ENABLE_OPENMP)
  find_package(OpenMP)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -DVIENNACL_WITH_OPENMP -DVIENNASHE_WITH_OPENMP")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DVIENNACL_WITH_OPENMP -DVIENNASHE_WITH_OPENMP")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  ADD_DEFINITIONS( -DVIENNACL_WITH_OPENMP -DVIENNASHE_WITH_OPENMP )
endif(ENABLE_OPENMP)

find_package(PETSc)


# ************************** Section 3: Miscellaneous compiler settings **************************

# Specify include directories
INCLUDE_DIRECTORIES(".")
INCLUDE_DIRECTORIES("external/")


# Set build type
IF(DEFINED CMAKE_BUILD_TYPE)
 SET (CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}")
ELSE()
 SET (CMAKE_BUILD_TYPE "Release")
ENDIF()


#
# Set high warning level on GCC and enable optimizations (or debug)
#
IF(CMAKE_COMPILER_IS_GNUCXX)
   IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
     ADD_DEFINITIONS( -Wall -O0 -g )           #debug build
   ELSE()
     ADD_DEFINITIONS( -Wall -pedantic -O3 )    #release build
   ENDIF()
ENDIF(CMAKE_COMPILER_IS_GNUCXX)


#
# Mac OS X specific linker part
#
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  IF(ENABLE_OPENCL)
    SET(CMAKE_EXE_LINKER_FLAGS " ${CMAKE_EXE_LINKER_FLAGS} -framework OpenCL")
  ENDIF(ENABLE_OPENCL)
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")


include_directories(
   ${PROJECT_BINARY_DIR}
   ${PROJECT_SOURCE_DIR}
   ${OPENCL_INCLUDE_DIRS}
   ${PETSC_INCLUDE_DIR})

# Set high warning level on GCC
if(ENABLE_PEDANTIC_FLAGS)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic")
endif()

IF(MSVC) #Visual Studio asks for the '/bigobj' flag
  ADD_DEFINITIONS(/bigobj)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
  SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /bigobj")
  SET(CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} /bigobj")
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /bigobj")
  SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} /bigobj")
ENDIF(MSVC)

### Turn logging on
IF (DISABLE_LOGGING)
  ADD_DEFINITIONS( -DVIENNASHE_LOG_DISABLE )
ENDIF (DISABLE_LOGGING)


IF(ENABLE_SERIAL_HDF5)
  INCLUDE_DIRECTORIES("external/hdf/include")

  ADD_DEFINITIONS( -DVIENNASHE_HAVE_HDF_SERIAL=1 )

  LINK_DIRECTORIES("${PROJECT_SOURCE_DIR}/external/hdf/lib")

  SET(HDF5_LIBRARIES "hdf5 -lz -lm" )

ENDIF(ENABLE_SERIAL_HDF5)


# Export
########

configure_file(cmake/FindOpenCL.cmake
   ${PROJECT_BINARY_DIR}/FindOpenCL.cmake COPYONLY)

configure_file(cmake/ViennaSHEConfig.cmake.in
   ${PROJECT_BINARY_DIR}/ViennaSHEConfig.cmake @ONLY)

configure_file(cmake/ViennaSHEConfigVersion.cmake.in
   ${PROJECT_BINARY_DIR}/ViennaSHEConfigVersion.cmake @ONLY)

if(CMAKE_MINOR_VERSION GREATER 6)  # export(PACKAGE ...) introduced with CMake 2.8.0
  export(PACKAGE ViennaSHE)
endif()

# Install
#########

install(FILES
   ${PROJECT_BINARY_DIR}/FindOpenCL.cmake
   ${PROJECT_BINARY_DIR}/ViennaSHEConfig.cmake
   ${PROJECT_BINARY_DIR}/ViennaSHEConfigVersion.cmake
   DESTINATION ${INSTALL_CMAKE_DIR} COMPONENT dev)
