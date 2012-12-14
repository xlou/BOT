CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(bot)

# Set path, etc.
SET(VIGRA_INSTALL_PATH "~/local" CACHE PATH "Additional install path for vigra")
SET(HDF5_INSTALL_PATH "~/local" CACHE PATH "Additional install path for HDF5")
SET(CPLEX_BASE_PATH "~/local" CACHE PATH "Additional install path for cplex")
SET(CPLEX_SYSTEM "x86-64_sles10_4.1" CACHE PATH "Specify the system of cplex")
SET(CPLEX_LIBFORMAT "static_pic" CACHE PATH "Specify the library format of cplex")

# Additional header files
INCLUDE_DIRECTORIES(
    ${VIGRA_INSTALL_PATH}/include
    ${HDF5_INSTALL_PATH}/include
    ${CPLEX_BASE_PATH}/cplex/include
    ${CPLEX_BASE_PATH}/concert/include
    ./include
    ./external
    ./external/object-features)

# Specify where to search for libraries
LINK_DIRECTORIES(
    ${VIGRA_INSTALL_PATH}/lib
    ${HDF5_INSTALL_PATH}/lib
    ${CPLEX_BASE_PATH}/cplex/lib/${CPLEX_SYSTEM}/${CPLEX_LIBFORMAT}
    ${CPLEX_BASE_PATH}/concert/lib/${CPLEX_SYSTEM}/${CPLEX_LIBFORMAT})

# Compiler flags
if(CMAKE_COMPILER_IS_GNUCXX)
    list(APPEND CMAKE_CXX_FLAGS -O3)
endif()
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -m64 -O -fPIC -fexceptions -DNDEBUG -DIL_STD")

FILE(GLOB BOT_ALL_SOURCE
    "./include/*.hxx"
    "./include/*.cxx")

ADD_LIBRARY(bot ${BOT_ALL_SOURCE})

# Add applications or tests
ADD_SUBDIRECTORY(apps)
SET(WITH_TESTS FALSE CACHE BOOL "Build test applications.")
if(WITH_TESTS)
    ADD_SUBDIRECTORY(tests)
endif()

# Installation
INSTALL(DIRECTORY include/ DESTINATION include/bot FILES_MATCHING PATTERN "*.hxx")
INSTALL(DIRECTORY external/object-features/ DESTINATION include/bot FILES_MATCHING PATTERN "*.hxx")
INSTALL(DIRECTORY external/SimpleIni/ DESTINATION include/bot/SimpleIni)
INSTALL(TARGETS bot DESTINATION lib)

