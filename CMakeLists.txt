cmake_minimum_required(VERSION 3.17)
project(alt_simulator)

set(CMAKE_CXX_STANDARD 17)

add_executable(alt_simulator main.cpp make_dataset.cpp alt_simulator.cpp)