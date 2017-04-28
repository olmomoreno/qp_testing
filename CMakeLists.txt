project(optimizer_all)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(src SRC_LIST)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -frounding-math")
include_directories(include)
link_directories(/home/olmo/Developer/snippets/optimizer_all/optimizer_all/lib)

add_executable(optimizer_all src/main.cpp)
target_link_libraries(optimizer_all cqp qpOASES CGAL gmp gsl gslcblas m)

