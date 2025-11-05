#include "tools.h"
#include <iostream>
#include "pch.h" 
#ifdef __linux__		
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif


Scalar mesh_scaling(Mesh& src_mesh, Mesh& tar_mesh)
{
    Vec3 max(-1e12, -1e12, -1e12);
    Vec3 min(1e12, 1e12, 1e12);
    for(auto it = src_mesh.vertices_begin(); it != src_mesh.vertices_end(); it++)
    {
        for(int j = 0; j < 3; j++)
        {
            if(src_mesh.point(*it)[j] > max[j])
            {
                max[j] = src_mesh.point(*it)[j];
            }
            if(src_mesh.point(*it)[j] < min[j])
            {
                min[j] = src_mesh.point(*it)[j];
            }
        }
    }

    for(auto it = tar_mesh.vertices_begin(); it != tar_mesh.vertices_end(); it++)
    {
        for(int j = 0; j < 3; j++)
        {
            if(tar_mesh.point(*it)[j] > max[j])
            {
                max[j] = tar_mesh.point(*it)[j];
            }
            if(tar_mesh.point(*it)[j] < min[j])
            {
                min[j] = tar_mesh.point(*it)[j];
            }
        }
    }
    Scalar scale = (max-min).norm();

    for(auto it = src_mesh.vertices_begin(); it != src_mesh.vertices_end(); it++)
    {
        Vec3 p = src_mesh.point(*it);
        p = p/scale;
        src_mesh.set_point(*it, p);
    }

    for(auto it = tar_mesh.vertices_begin(); it != tar_mesh.vertices_end(); it++)
    {
        Vec3 p = tar_mesh.point(*it);
        p = p/scale;
        tar_mesh.set_point(*it, p);
    }

    return scale;
}



Vector3 Vec2Eigen(Vec3 s)
{
    return Vector3(s[0], s[1], s[2]);
}


#ifdef __linux__
bool my_mkdir(std::string file_path)
{
    if(access(file_path.c_str(), 06))
   {
       std::cout << "file_path : (" << file_path << ") didn't exist or no write ability!!" << std::endl;
       if(mkdir(file_path.c_str(), S_IRWXU))
       {
           std::cout << "mkdir " << file_path << " is wrong! please check upper path " << std::endl;
           exit(0);
       }
       std::cout<< "mkdir " << file_path << " success!! " << std::endl;
   }
}
#endif
