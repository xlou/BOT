CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(bot-apps)

# test: bot tracking learner
ADD_EXECUTABLE(TrackingTrainer-hdf5 TrackingTrainer-hdf5.cxx ${BOT_ALL_SOURCE})
# the following "-lm -pthread" is needed for cplex
SET_TARGET_PROPERTIES(TrackingTrainer-hdf5 PROPERTIES LINK_FLAGS "-lm -pthread")
IF(USE_TIFF)
    TARGET_LINK_LIBRARIES(TrackingTrainer-hdf5 tiff vigraimpex ilocplex cplex concert)
ELSE(USE_TIFF)
    TARGET_LINK_LIBRARIES(TrackingTrainer-hdf5 hdf5 hdf5_hl vigraimpex ilocplex cplex concert)
ENDIF(USE_TIFF)

# test: bot tracking predictor
ADD_EXECUTABLE(TrackingPredictor-hdf5 TrackingPredictor-hdf5.cxx ${BOT_ALL_SOURCE})
# the following "-lm -pthread" is needed for cplex
SET_TARGET_PROPERTIES(TrackingPredictor-hdf5 PROPERTIES LINK_FLAGS "-lm -pthread") 
IF(USE_TIFF)
    TARGET_LINK_LIBRARIES(TrackingPredictor-hdf5 tiff vigraimpex ilocplex cplex concert)
ELSE(USE_TIFF)
    TARGET_LINK_LIBRARIES(TrackingPredictor-hdf5 hdf5 hdf5_hl vigraimpex ilocplex cplex concert)
ENDIF(USE_TIFF)
