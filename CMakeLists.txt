CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(bot)

# Use hdf5 or tiff
OPTION(USE_TIFF "Set to ON to use tiff as input data format." ON)
IF(USE_TIFF)
    ADD_DEFINITIONS(-DUSE_TIFF)
ENDIF(USE_TIFF)

OPTION(USE_HDF5 "Set to ON to use hdf5 as input data format." OFF)
IF(USE_HDF5)
    ADD_DEFINITIONS(-DUSE_HDF5)
ENDIF(USE_HDF5)

# Set path, etc.
SET(VIGRA_INSTALL_PATH "" CACHE PATH "Additional install path for vigra")
SET(LIBTIFF_INSTALL_PATH "" CACHE PATH "Libtiff installation path")
SET(HDF5_INSTALL_PATH "" CACHE PATH "Additional install path for HDF5")
SET(CPLEX_BASE_PATH "" CACHE PATH "Additional install path for cplex")
SET(CPLEX_SYSTEM "x86-64_sles10_4.1" CACHE PATH "Specify the system of cplex")
SET(CPLEX_LIBFORMAT "static_pic" CACHE PATH "Specify the library format of cplex")

# Additional header files
INCLUDE_DIRECTORIES(
    ${VIGRA_INSTALL_PATH}/include
    ${LIBTIFF_INSTALL_PATH}/include
    ${HDF5_INSTALL_PATH}/include
    ${CPLEX_BASE_PATH}/cplex/include
    ${CPLEX_BASE_PATH}/concert/include
    ./include
    ./external
    ./external/object-features)

# Specify where to search for libraries
LINK_DIRECTORIES(
    ${VIGRA_INSTALL_PATH}/lib
    ${LIBTIFF_INSTALL_PATH}/lib
    ${HDF5_INSTALL_PATH}/lib
    ${CPLEX_BASE_PATH}/cplex/lib/${CPLEX_SYSTEM}/${CPLEX_LIBFORMAT}
    ${CPLEX_BASE_PATH}/concert/lib/${CPLEX_SYSTEM}/${CPLEX_LIBFORMAT})

# Sometime library names are different between Win and Linux
IF(WIN32)
        SET(LIB_PREFIX "lib")
ELSE(WIN32)
        SET(LIB_PREFIX "")
ENDIF(WIN32)

# Compiler flags
IF(CMAKE_COMPILER_IS_GNUCXX)
    list(APPEND CMAKE_CXX_FLAGS -O3)
ENDIF()
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -m64 -O -fPIC -fexceptions -DNDEBUG -DIL_STD")

FILE(GLOB BOT_ALL_SOURCE
    "./include/*.hxx"
    "./include/*.cxx")

ADD_LIBRARY(bot ${BOT_ALL_SOURCE})

# Add applications or tests
ADD_SUBDIRECTORY(apps)
SET(WITH_TESTS FALSE CACHE BOOL "Build test applications.")
IF(WITH_TESTS)
    ADD_SUBDIRECTORY(tests)
ENDIF()

# Installation
INSTALL(DIRECTORY include/ DESTINATION include/bot FILES_MATCHING PATTERN "*.hxx")
INSTALL(DIRECTORY external/object-features/ DESTINATION include/bot FILES_MATCHING PATTERN "*.hxx")
INSTALL(DIRECTORY external/SimpleIni/ DESTINATION include/bot/SimpleIni)
INSTALL(TARGETS bot DESTINATION lib)

