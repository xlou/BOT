CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(bot-tests)

#OPTION(BUILD_SHARED_LIBS OFF)

SET(VIGRA_INSTALL_PATH "~/local" CACHE PATH "Additional install path for vigra")
SET(HDF5_INSTALL_PATH "~/local" CACHE PATH "Additional install path for HDF5")
SET(VIGRANUMPY_EXTENSIONS_INSTALL_PATH "~/local" CACHE PATH "Additional install path for vigranumpy-extensions")
SET(CPLEX_BASE_PATH "~/local" CACHE PATH "Additional install path for cplex")
SET(CPLEX_SYSTEM "x86-64_sles10_4.1" CACHE PATH "Specify the system of cplex")
SET(CPLEX_LIBFORMAT "static_pic" CACHE PATH "Specify the library format of cplex")

INCLUDE_DIRECTORIES(
    ${VIGRA_INSTALL_PATH}/include
    ${HDF5_INSTALL_PATH}/include
    ${VIGRANUMPY_EXTENSIONS_INSTALL_PATH}
    ${CPLEX_BASE_PATH}/cplex/include
    ${CPLEX_BASE_PATH}/concert/include
    ../include
    ../external)

# Specify where to search for libraries
LINK_DIRECTORIES(
    ${VIGRA_INSTALL_PATH}/lib
    ${HDF5_INSTALL_PATH}/lib
    ${CPLEX_BASE_PATH}/cplex/lib/${CPLEX_SYSTEM}/${CPLEX_LIBFORMAT}
    ${CPLEX_BASE_PATH}/concert/lib/${CPLEX_SYSTEM}/${CPLEX_LIBFORMAT})

file(GLOB BOT_ALL_SOURCE
    "../include/*.hxx"
    "../include/*.cxx")

# vigranumpy-extensions library
ADD_LIBRARY(vigranumpy-extensions
    ${VIGRANUMPY_EXTENSIONS_INSTALL_PATH}/objectFeatures.cxx)

# test: c++ factory method
ADD_EXECUTABLE(test-cpp-factory
    test-cpp-factory.cxx
    Stooge.hxx
    Moe.hxx
    Larry.hxx
    Curly.hxx
    Stooge.cxx)

# test: vigranumpy-extensions features
ADD_EXECUTABLE(test-vigranumpy-extensions
    test-vigranumpy-extensions.cxx)
#TARGET_LINK_LIBRARIES(test-vigranumpy-extensions vigranumpy-extensions vigraimpex)

# test: bot tracking learner
ADD_EXECUTABLE(test-TrackingTrainer
    test-TrackingTrainer.cxx
    ${BOT_ALL_SOURCE})
SET_TARGET_PROPERTIES(test-TrackingTrainer PROPERTIES LINK_FLAGS "-lm -pthread")
TARGET_LINK_LIBRARIES(test-TrackingTrainer hdf5 hdf5_hl vigraimpex ilocplex cplex concert)

# test: bot tracking predictor
ADD_EXECUTABLE(test-TrackingPredictor
    test-TrackingPredictor.cxx
    ${BOT_ALL_SOURCE})
SET_TARGET_PROPERTIES(test-TrackingPredictor PROPERTIES LINK_FLAGS "-lm -pthread")
TARGET_LINK_LIBRARIES(test-TrackingPredictor hdf5 hdf5_hl vigraimpex ilocplex cplex concert)

# test: bot hypothesis space
ADD_EXECUTABLE(test-HypothesisSpace
    test-HypothesisSpace.cxx
    ${BOT_ALL_SOURCE})
TARGET_LINK_LIBRARIES(test-HypothesisSpace hdf5 hdf5_hl vigraimpex)

# test: bot event configuration
ADD_EXECUTABLE(test-EventConfiguration
    test-EventConfiguration.cxx
    ${BOT_ALL_SOURCE})

# test: bot object/party list generator
ADD_EXECUTABLE(test-ObjectsGenerator
    test-ObjectsGenerator.cxx
    ${BOT_ALL_SOURCE})
TARGET_LINK_LIBRARIES(test-ObjectsGenerator hdf5 hdf5_hl vigraimpex)

# test: bot nearest neighbor generator
ADD_EXECUTABLE(test-NearestNeighborGenerator
    test-NearestNeighborGenerator.cxx
    ${BOT_ALL_SOURCE})

# test: bot nearest neighbor generator
ADD_EXECUTABLE(test-ConvexImage
    test-ConvexImage.cxx
    ${BOT_ALL_SOURCE})
    
# test: bot object feature factory
#ADD_EXECUTABLE(test-ObjectFeatureFactory
#    test-ObjectFeatureFactory.cxx
#    ${BOT_ALL_SOURCE})

# Compiler flags
if(CMAKE_COMPILER_IS_GNUCXX)
    list(APPEND CMAKE_CXX_FLAGS -O3)
endif()

# test: cplex ilp/qp solver
ADD_EXECUTABLE(test-CplexSolver
    test-CplexSolver.cxx)
SET_TARGET_PROPERTIES(test-CplexSolver PROPERTIES LINK_FLAGS "-lm -pthread")
TARGET_LINK_LIBRARIES(test-CplexSolver ilocplex cplex concert)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -m64 -O -fPIC -fexceptions -DNDEBUG -DIL_STD")