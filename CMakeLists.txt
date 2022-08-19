cmake_minimum_required(VERSION 3.5)
project(ab-mrp-cutplan)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_STD -D_GLIBCXX_USE_CXX11_ABI=0 -std=c++0x ")

include_directories(/opt/gurobi952/linux64/include /opt/boost_1_79_0 /opt/lemon-1.3.1/build /opt/lemon-1.3.1)

link_directories(/opt/gurobi952/linux64/lib /opt/boost_1_79_0 /opt/lemon-1.3.1/build/lemon)

set(SOURCE_FILES main.cpp headers/Graph.h src/Graph.cpp headers/Include.h src/Model.cpp headers/Model.h src/Arc.cpp headers/Arc.h)

add_executable(abmodel ${SOURCE_FILES})

target_link_libraries(abmodel gurobi_c++ gurobi95 pthread m emon)
