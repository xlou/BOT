cmake_minimum_required(VERSION 2.6)
project(bot)

include_directories(external/ external/object-features)

file(GLOB bot_srcs include/*.cxx)
add_library(bot ${bot_srcs})

set(WITH_TESTS TRUE CACHE BOOL "Build test applications.")
if(WITH_TESTS)
    add_subdirectory(tests)
endif()

##
# installation
##
# header files
install(DIRECTORY include/ DESTINATION include/bot FILES_MATCHING PATTERN "*.hxx")
install(DIRECTORY external/object-features/ DESTINATION include/bot FILES_MATCHING PATTERN "*.hxx")
install(DIRECTORY external/SimpleIni/ DESTINATION include/bot/SimpleIni)

# lib
install(TARGETS bot DESTINATION lib)

