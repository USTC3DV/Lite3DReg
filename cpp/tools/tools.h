#ifndef TOOL_H_
#define TOOL_H_
#include "Types.h"
// Convert Mesh to libigl format to calculate geodesic distance

Scalar mesh_scaling(Mesh& src_mesh, Mesh& tar_mesh);
Vector3 Vec2Eigen(Vec3 s);
#ifdef __linux__
bool my_mkdir(std::string file_path);
#endif

#endif
