CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(bot-tests)

# test: c++ factory method
ADD_EXECUTABLE(test-cpp-factory
    test-cpp-factory.cxx
    Stooge.hxx
    Moe.hxx
    Larry.hxx
    Curly.hxx
    Stooge.cxx)

# test tiff file reading
IF(USE_TIFF)
    ADD_EXECUTABLE(test-TIFFReaderWriter test-TIFFReaderWriter.cxx ${BOT_ALL_SOURCE})
    TARGET_LINK_LIBRARIES(test-TIFFReaderWriter ${LIB_PREFIX}tiff)
ENDIF(USE_TIFF)

# test: bot hypothesis space
IF(USE_HDF5)
    ADD_EXECUTABLE(test-HypothesisSpace test-HypothesisSpace.cxx ${BOT_ALL_SOURCE})
    TARGET_LINK_LIBRARIES(test-HypothesisSpace hdf5 hdf5_hl vigraimpex)
ENDIF(USE_HDF5)

# test: bot event configuration
ADD_EXECUTABLE(test-EventConfiguration
    test-EventConfiguration.cxx
    ${BOT_ALL_SOURCE})

# test: bot object/party list generator
IF(USE_HDF5)
    ADD_EXECUTABLE(test-ObjectsGenerator test-ObjectsGenerator.cxx ${BOT_ALL_SOURCE})
    TARGET_LINK_LIBRARIES(test-ObjectsGenerator hdf5 hdf5_hl vigraimpex)
ENDIF(USE_HDF5)

# test: bot nearest neighbor generator
ADD_EXECUTABLE(test-NearestNeighborGenerator test-NearestNeighborGenerator.cxx ${BOT_ALL_SOURCE})

# test: bot nearest neighbor generator
ADD_EXECUTABLE(test-ConvexImage test-ConvexImage.cxx ${BOT_ALL_SOURCE})
    
# test: cplex ilp/qp solver
ADD_EXECUTABLE(test-CplexSolver
    test-CplexSolver.cxx)
SET_TARGET_PROPERTIES(test-CplexSolver PROPERTIES LINK_FLAGS "-lm -pthread")
TARGET_LINK_LIBRARIES(test-CplexSolver tiff ilocplex cplex concert)

