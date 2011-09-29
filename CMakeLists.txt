cmake_minimum_required(VERSION 2.6)
project(bot)

##################
### components ###
##################
set(WITH_TESTS TRUE CACHE BOOL "Build test applications.")

if(WITH_TESTS)
    add_subdirectory(tests)
endif()

install(DIRECTORY include/ DESTINATION include/bot FILES_MATCHING PATTERN "*.hxx")