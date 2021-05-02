
PROJECT(ViennaSHE)

cmake_minimum_required(VERSION 2.6)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

set(VERSION_MAJOR 1)
set(VERSION_MINOR 2)
set(VERSION_PATCH 0)
set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})
<<<<<<< Updated upstream


# If you want the full compiler output, uncomment the following line
#SET(CMAKE_VERBOSE_MAKEFILE ON)


# ************************** Section 1: Configure  **************************

#
# Uncomment the following line to enable the GPU-accelerated OpenCL-based solver:
#
#SET(ENABLE_OPENCL ON)


#
###################### END OF USER CUSTOMIZATION OPTIONS ################################
#
# In most situations, there is no need to change anything below this line.
#
=======


# If you want the full compiler output, uncomment the following line
#SET(CMAKE_VERBOSE_MAKEFILE ON)


# ************************** Section 1: Configure  **************************

#
# Uncomment the following line to enable the GPU-accelerated OpenCL-based solver:
#
#SET(ENABLE_OPENCL ON)

>>>>>>> Stashed changes

#
###################### END OF USER CUSTOMIZATION OPTIONS ################################
#
# In most situations, there is no need to change anything below this line.
#

<<<<<<< Updated upstream
# **************************  Section 2: Configure CMake **************************


list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include(ViennaSHECommon)

# Set release build:
IF(DEFINED CMAKE_BUILD_TYPE)
 SET (CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE})
ELSE()
 SET (CMAKE_BUILD_TYPE Release)
ENDIF()


include_directories(BEFORE ${PROJECT_SOURCE_DIR})
=======

# **************************  Section 2: Configure CMake **************************


list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include(ViennaSHECommon)

# Set release build:
IF(DEFINED CMAKE_BUILD_TYPE)
 SET (CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE})
ELSE()
 SET (CMAKE_BUILD_TYPE Release)
ENDIF()

>>>>>>> Stashed changes

add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(python)

if(BUILD_TESTING)
   INCLUDE(CTest)
   add_subdirectory(tests)
endif()


add_subdirectory(libviennashe)


# Install
#########

install(DIRECTORY viennashe
   DESTINATION ${INSTALL_INCLUDE_DIR} COMPONENT dev
   FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")

