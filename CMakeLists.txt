cmake_minimum_required(VERSION 3.16)
project(Lightweight3DRegLib)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# 第三方库路径
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# --- 第三方库
find_package(Eigen3 3.4.0 REQUIRED)
message(STATUS "Found Eigen3 version: ${EIGEN3_VERSION_STRING}")
find_package(OpenMesh REQUIRED)
find_package(OpenMP REQUIRED)

# --- 添加 pybind11 支持 ---
find_package(pybind11 CONFIG REQUIRED)

option(USEPARDISO "Use Pardiso solver" ON)

# --- 源文件
file(GLOB_RECURSE SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/cpp/core/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/cpp/io/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/cpp/tools/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/cpp/tools/geodesic/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/cpp/tools/preprocess/*.cpp"
)

file(GLOB_RECURSE HEADERS
    "${CMAKE_CURRENT_SOURCE_DIR}/cpp/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/cpp/core/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/cpp/io/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/cpp/tools/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/cpp/tools/geodesic/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/cpp/tools/preprocess/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/*.h"
)

# --- 创建静态库
add_library(Lightweight3DRegLib STATIC ${SOURCES} ${HEADERS})

target_precompile_headers(Lightweight3DRegLib PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/cpp/pch.h")

set_target_properties(Lightweight3DRegLib PROPERTIES POSITION_INDEPENDENT_CODE ON)

# --- include 路径
target_include_directories(Lightweight3DRegLib
    PUBLIC 
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/core"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/io"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/tools"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/tools/geodesic"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/tools/preprocess"
        "${CMAKE_CURRENT_SOURCE_DIR}/3rdparty"
        ${EIGEN3_INCLUDE_DIRS}
)

# --- 链接库
target_link_libraries(Lightweight3DRegLib
    PUBLIC OpenMeshCore OpenMeshTools
    PUBLIC OpenMP::OpenMP_CXX
)

# --- Pardiso / MKL
if(USEPARDISO AND ${CMAKE_HOST_SYSTEM_NAME} MATCHES "Linux")
    if(EXISTS /opt/intel/oneapi/mkl/latest/include/mkl_pardiso.h)
        target_compile_definitions(Lightweight3DRegLib PUBLIC USE_PARDISO)
        target_include_directories(Lightweight3DRegLib SYSTEM PUBLIC "/opt/intel/oneapi/mkl/latest/include")
        target_compile_options(Lightweight3DRegLib PUBLIC "-m64")
        target_link_libraries(Lightweight3DRegLib PRIVATE
            "-Wl,--start-group /opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_intel_lp64.a"
            "/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_gnu_thread.a"
            "/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_core.a"
            "-Wl,--end-group -lgomp -lpthread -lm -ldl"
        )
    else()
        message("Cannot find Intel MKL at /opt/intel/oneapi/mkl. Using Eigen LDLT.")
    endif()
endif()

# --- 添加 Python 模块 ---
pybind11_add_module(pyregister ${CMAKE_CURRENT_SOURCE_DIR}/python/pyregister.cpp)

target_link_libraries(pyregister 
    PRIVATE 
        Lightweight3DRegLib  
        OpenMeshCore 
        OpenMeshTools
        OpenMP::OpenMP_CXX
)

target_include_directories(pyregister 
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/core"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/io"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/tools"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/tools/geodesic"
        "${CMAKE_CURRENT_SOURCE_DIR}/cpp/tools/preprocess"
        "${CMAKE_CURRENT_SOURCE_DIR}/3rdparty"
        ${EIGEN3_INCLUDE_DIRS}
)

target_precompile_headers(pyregister REUSE_FROM Lightweight3DRegLib)

if(USEPARDISO AND ${CMAKE_HOST_SYSTEM_NAME} MATCHES "Linux")
    if(EXISTS /opt/intel/oneapi/mkl/latest/include/mkl_pardiso.h)
        target_compile_definitions(pyregister PRIVATE USE_PARDISO)
        target_include_directories(pyregister SYSTEM PRIVATE "/opt/intel/oneapi/mkl/latest/include")
        target_compile_options(pyregister PRIVATE "-m64")
        target_link_libraries(pyregister PRIVATE
            "-Wl,--start-group /opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_intel_lp64.a"
            "/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_gnu_thread.a"
            "/opt/intel/oneapi/mkl/latest/lib/intel64/libmkl_core.a"
            "-Wl,--end-group -lgomp -lpthread -lm -ldl"
        )
    endif()
endif()

set_target_properties(pyregister PROPERTIES 
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python
)


# =====================================================================
# --- 添加示例可执行程序 (examples/c++/main.cpp) ---
# =====================================================================
file(GLOB EXAMPLE_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/examples/c++/*.cpp"
)


add_executable(Lightweight3DRegDemo ${EXAMPLE_SOURCES})

target_precompile_headers(Lightweight3DRegDemo PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/cpp/pch.h")


target_link_libraries(Lightweight3DRegDemo
    PRIVATE
        Lightweight3DRegLib
        OpenMeshCore 
        OpenMeshTools
        OpenMP::OpenMP_CXX
)



set_target_properties(Lightweight3DRegDemo PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/examples
)