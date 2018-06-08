cmake_minimum_required(VERSION 3.10)
project(Valentin)

set(CMAKE_CXX_STANDARD 11)


find_package(Boost COMPONENTS system filesystem REQUIRED)



#include_directories(/home/lamure/Documents/developpement/boost/)

SET(CMAKE_C_COMPILER g++)
SET(CMAKE_CXX_COMPILER g++)
set(SOURCE_FILES main.txt main2.txt mainv2.cpp emd_test.cpp emd_test.hpp)
add_executable(Valentin ${SOURCE_FILES})
link_directories(/usr/lib64/)


target_link_libraries(Valentin boost_graph boost_system)
target_link_libraries(Valentin boost_thread-mt)