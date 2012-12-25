CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(bot-apps)

# test: bot tracking learner
ADD_EXECUTABLE(TrackingTrainer TrackingTrainer.cxx ${BOT_ALL_SOURCE})
# the following "-lm -pthread" is needed for cplex
SET_TARGET_PROPERTIES(TrackingTrainer PROPERTIES LINK_FLAGS "-lm -pthread")
IF(USE_TIFF)
    TARGET_LINK_LIBRARIES(TrackingTrainer ${LIB_PREFIX}tiff vigraimpex ilocplex concert ${CPLEX_LIBNAME})
ELSE(USE_TIFF)
    TARGET_LINK_LIBRARIES(TrackingTrainer hdf5 hdf5_hl vigraimpex ilocplex concert ${CPLEX_LIBNAME})
ENDIF(USE_TIFF)

# test: bot tracking predictor
ADD_EXECUTABLE(TrackingPredictor TrackingPredictor.cxx ${BOT_ALL_SOURCE})
# the following "-lm -pthread" is needed for cplex
SET_TARGET_PROPERTIES(TrackingPredictor PROPERTIES LINK_FLAGS "-lm -pthread") 
IF(USE_TIFF)
    TARGET_LINK_LIBRARIES(TrackingPredictor ${LIB_PREFIX}tiff vigraimpex ilocplex concert ${CPLEX_LIBNAME})
ELSE(USE_TIFF)
    TARGET_LINK_LIBRARIES(TrackingPredictor hdf5 hdf5_hl vigraimpex ilocplex concert ${CPLEX_LIBNAME})
ENDIF(USE_TIFF)
