# - Try to find OpenMesh
# This is a minimal FindOpenMesh.cmake for systems with libopenmesh-dev installed

# 查找头文件
find_path(OPENMESH_INCLUDE_DIR OpenMesh/Core/IO/MeshIO.hh
          PATHS /usr/include /usr/include/OpenMesh
          PATH_SUFFIXES OpenMesh)

# 查找库文件
find_library(OPENMESH_CORE_LIBRARY OpenMeshCore
             PATHS /usr/lib /usr/lib/x86_64-linux-gnu)

find_library(OPENMESH_TOOLS_LIBRARY OpenMeshTools
             PATHS /usr/lib /usr/lib/x86_64-linux-gnu)

# 设置变量给 CMake 使用
set(OPENMESH_LIBRARIES ${OPENMESH_CORE_LIBRARY} ${OPENMESH_TOOLS_LIBRARY})
set(OPENMESH_FOUND TRUE)

# 输出信息
if(OPENMESH_FOUND)
    message(STATUS "Found OpenMesh:")
    message(STATUS "  Include dir: ${OPENMESH_INCLUDE_DIR}")
    message(STATUS "  Libraries: ${OPENMESH_LIBRARIES}")
else()
    message(FATAL_ERROR "OpenMesh not found")
endif()
