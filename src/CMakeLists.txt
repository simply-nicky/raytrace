file(GLOB SOURCES *.cpp)
add_library(sources STATIC ${SOURCES})
target_link_libraries(sources ${Boost_LIBRARIES})