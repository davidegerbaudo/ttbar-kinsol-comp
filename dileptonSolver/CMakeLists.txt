cmake_minimum_required(VERSION 2.6)

project(test_solver) 

#options here
set(CMAKE_CXX_FLAGS "-Wall -march=native")

#get ROOT stuff here
exec_program(root-config ARGS "--incdir" OUTPUT_VARIABLE ROOT_INCLUDE_DIR )
include_directories(include ${ROOT_INCLUDE_DIR})
#exec_program(root-config ARGS "--libdir" OUTPUT_VARIABLE ROOT_LIBRARY_DIR )
#link_directories(${ROOT_LIBRARY_DIR})
exec_program(root-config ARGS "--libs" OUTPUT_VARIABLE ROOT_LIBS )
#message(STATUS ${ROOT_LIBS})

#project here
set(test_solver_SRCS example.cc DileptonAnalyticalSolver.cc)
add_executable(test_solver ${test_solver_SRCS})
