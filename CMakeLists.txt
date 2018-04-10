cmake_minimum_required(VERSION 3.6 FATAL_ERROR)

project(raytrace_build VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
# set(Boost_DEBUG 1)


# if(NOT WIN32)
# file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/Al2O3.txt 
# DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
# endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if (NOT DEFINED     BOOST_ROOT        AND
NOT DEFINED ENV{BOOST_ROOT}       AND
NOT DEFINED     BOOST_LIBRARYDIR  AND
NOT DEFINED ENV{BOOST_LIBRARYDIR})
if (APPLE)
set(BOOST_ROOT "/usr/local/boost_1_66_0")
set(BOOST_LIBRARYDIR "/usr/local/lib")
elseif (WIN32)
set(BOOST_ROOT "C:/'Program Files (x86)'/boost/boost_1_66_0")
set(BOOST_LIBRARYDIR "C:/'Program Files (x86)'/boost/boost_1_66_0/stage/lib")
elseif(CYGWIN)
set(BOOST_ROOT "C:\\cygwin64\\usr\\local\\include")
set(BOOST_LIBRARYDIR "C:\\cygwin64\\usr\\local\\lib")
elseif(UNIX AND NOT APPLE)
set(BOOST_ROOT "/usr/local/include")
set(BOOST_LIBRARYDIR "/usr/local/lib")
endif()
endif()
if (UNIX OR WIN32)
set(Boost_USE_STATIC_LIBS TRUE)
endif()
find_package(Boost 1.66 REQUIRED COMPONENTS filesystem program_options system)

if(Boost_FOUND)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
add_subdirectory(src)
add_executable(raytrace main.cpp)
target_link_libraries(raytrace ${Boost_LIBRARIES} sources)
add_executable(Test test.cpp)
target_link_libraries(Test sources)
add_custom_command(
                TARGET raytrace POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy
                        ${CMAKE_SOURCE_DIR}/Al2O3.txt 
                        ${CMAKE_CURRENT_BINARY_DIR}/Al2O3.txt
                 )
if(UNIX AND NOT APPLE)
find_package(Threads REQUIRED)
target_link_libraries(raytrace Threads::Threads)
target_link_libraries(Test Threads::Threads)
endif()
endif()

