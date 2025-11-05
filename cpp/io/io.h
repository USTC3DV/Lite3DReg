#ifndef IO_H
#define IO_H
#include <cstdio>
#include "pch.h"
#include "Types.h"


inline bool read_by_openmesh(const std::string filename, Mesh& mesh)
{
    OpenMesh::IO::Options opt_read = OpenMesh::IO::Options::VertexNormal;
    mesh.request_vertex_normals();
    bool read_OK = OpenMesh::IO::read_mesh(mesh, filename,opt_read);

	std::cout << "filename = " << filename << std::endl;
    if (read_OK)
    {
        mesh.request_vertex_status();
        mesh.request_edge_status();
        mesh.request_face_status();

        mesh.request_face_normals();

        mesh.update_face_normals();
        if(mesh.n_faces()>0)
            mesh.update_vertex_normals();

        Vec3 MeshScales;
        MeshScales[0] = 0.0; MeshScales[1] = 0.0; MeshScales[2] = 0.0;
        for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
        {
            MeshScales += mesh.point(*v_it);
        }
        MeshScales /= mesh.n_vertices();
        return true;
    }
    std::cout << "#vertices = " << mesh.n_vertices() << std::endl;
    return false;
}

//---openmesh for output 
inline bool write_by_openmesh(const char* filename, Mesh& mesh, Scalar scale)
{
    for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
    {

       mesh.point(*v_it) = mesh.point(*v_it)*scale;
    }
    OpenMesh::IO::Options opt_write = OpenMesh::IO::Options::VertexNormal;
    bool ok = OpenMesh::IO::write_mesh(mesh, filename,opt_write);
	return ok;
}

// If it exists, read the initial rigid transformation matrix.
template <class MatrixType>
bool read_transMat(MatrixType& trans, const std::string& filename)
{
	std::ifstream input(filename);
	std::string line;
	int rows, cols;
	std::vector<std::vector<double>> total_data;
	while (getline(input, line)) {
        if(line[0] == 'V' || line[0] == 'M')
            continue;
		std::istringstream iss(line);
		std::vector<double> lineVec;
		while (iss) {
			double item;
			if (iss >> item)
				lineVec.push_back(item);
		}
		cols = lineVec.size();
		total_data.push_back(lineVec);
	}
	if (total_data.size() == 0)
	{
		std::cout << filename << " is empty !! " << std::endl;
		return false;
	}
	rows = total_data.size();
	trans.resize(rows, cols);
    std::cout << "rows = " << rows << " cols = " << cols << std::endl;
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			trans(i, j) = total_data[i][j];
		}
	}
	input.close();
    std::cout << "read trans = \n" << trans << std::endl;
	return true;
}

#endif // IO_H