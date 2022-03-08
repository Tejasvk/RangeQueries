cmake_minimum_required(VERSION 2.8.7)
project(rangequeries)
set(CMAKE_BUILD_TYPE Release)

#Specify C++11 flag
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX14)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O3")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -O3")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++14 support. Please use a different C++ compiler.")
endif()


#However, the file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "sources/*.cpp")
file(GLOB INCLUDES "Headers/*.h")
#Generate the shared library from the sources
add_library(rangequeries STATIC ${SOURCES})

#Set the location for library installation -- i.e., /usr/lib in this case
# not really necessary in this example. Use "sudo make install" to apply
#install(TARGETS rangequeries DESTINATION /usr/lib)
#install(FILES src/rangequeries DESTINATION /usr/local/include)
#install(FILES ${INCLUDES} DESTINATION /usr/local/include/streamingcc_include)



