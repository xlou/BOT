CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(bot-apps)

# test: bot tracking learner
ADD_EXECUTABLE(TrackingTrainer
    TrackingTrainer.cxx
    ${BOT_ALL_SOURCE})
SET_TARGET_PROPERTIES(TrackingTrainer PROPERTIES LINK_FLAGS "-lm -pthread")
TARGET_LINK_LIBRARIES(TrackingTrainer hdf5 hdf5_hl vigraimpex ilocplex cplex concert)

# test: bot tracking predictor
ADD_EXECUTABLE(TrackingPredictor
    TrackingPredictor.cxx
    ${BOT_ALL_SOURCE})
SET_TARGET_PROPERTIES(TrackingPredictor PROPERTIES LINK_FLAGS "-lm -pthread")
TARGET_LINK_LIBRARIES(TrackingPredictor hdf5 hdf5_hl vigraimpex ilocplex cplex concert)

