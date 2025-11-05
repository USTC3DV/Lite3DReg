cmake_minimum_required(VERSION 3.16)
project(MyExternalApp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()




add_subdirectory(002hugging_face)
add_executable(MyApp main.cpp)
target_link_libraries(MyApp PRIVATE Lightweight3DRegLib)

