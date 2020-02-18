cmake_minimum_required(VERSION 3.14)
project(GeneticAlgorithm C)

set(CMAKE_C_STANDARD 99)

add_executable(GeneticAlgorithm main.c RKF78.h RKF78.c fitness.c fitness.h GACore.h GACore.c)
target_link_libraries(GeneticAlgorithm m)