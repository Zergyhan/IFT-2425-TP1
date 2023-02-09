cmake_minimum_required(VERSION 3.11)
project(TP1 C)

set(CMAKE_C_STANDARD 17)
# Add executable called "tp1" that is built from the source files and linked with the libraries
add_executable(tp1 main.c)
target_link_libraries(tp1 m)
