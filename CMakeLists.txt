cmake_minimum_required(VERSION 3.10)
project(Valentin)

set(CMAKE_CXX_STANDARD 11)


find_package(Boost COMPONENTS system filesystem REQUIRED)



#include_directories(/home/lamure/Documents/developpement/boost/)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


SET(CMAKE_C_COMPILER g++)
SET(CMAKE_CXX_COMPILER g++)
set(SOURCE_FILES main.txt main2.txt mainv2.cpp emd_test.cpp emd.h)
add_executable(Valentin ${SOURCE_FILES})
link_directories(/usr/lib64/)


target_link_libraries(Valentin boost_graph boost_system)
target_link_libraries(Valentin boost_thread-mt)