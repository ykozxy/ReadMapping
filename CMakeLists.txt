cmake_minimum_required(VERSION 3.25)

project(ReadMapping)

set(CMAKE_CXX_STANDARD 17)

file(GLOB SOURCES "src/*.cpp")
file(GLOB HEADERS "src/*.h")
#file(GLOB INCLUDES "include/*.cpp" "include/*.h")

add_executable(debug ${SOURCES})
target_include_directories(debug PRIVATE include)
target_compile_options(debug PRIVATE -Wall -Wextra -g -O0)

add_executable(release ${SOURCES})
target_include_directories(release PRIVATE include)
target_compile_options(release PRIVATE -Wall -Wextra -O3)
