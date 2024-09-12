cmake_minimum_required(VERSION 3.10)
project(hello_world)
 
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    /usr/local/include
    /usr/include/
    /usr/include/eigen3
)
 
set(CMAKE_BUILD_TYPE "Debug") # 默认是Release模式，设置为Debug才能调试
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin) # 设置可执行文件的输出目录
 
add_executable(hello_world test_eigen.cpp)