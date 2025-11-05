import os
import numpy as np
import trimesh
from plyfile import PlyData
import pyregister


def read_vertices_normals(file_path):
    """
    Read vertices and normals from a point cloud or mesh file.
    Supports formats such as .ply, .obj, .xyz, etc.
    """
    ext = os.path.splitext(file_path)[1].lower()

    if ext == ".ply":
        plydata = PlyData.read(file_path)
        vertex_data = plydata['vertex']
        vertices = np.vstack([vertex_data['x'], vertex_data['y'], vertex_data['z']]).astype(np.float64)

        normals = None
        if {'nx', 'ny', 'nz'}.issubset(vertex_data.data.dtype.names):
            normals = np.vstack([vertex_data['nx'], vertex_data['ny'], vertex_data['nz']]).astype(np.float64)

    else:
        mesh = trimesh.load(file_path, process=False)
        vertices = mesh.vertices.T.astype(np.float64)
        normals = getattr(mesh, "vertex_normals", None)
        if normals is not None and len(normals) > 0:
            normals = normals.T.astype(np.float64)
        else:
            normals = np.empty((0, 0))

    return vertices, normals


def run_rigid_file(file_target, file_source, output_path):
    """
    Perform rigid registration (FRICP method) by directly reading files.

    Parameters
    ----------
    file_target : str
        Path to the target point cloud or mesh file.
    file_source : str
        Path to the source point cloud or mesh file.
    output_path : str
        Output file path for the registered result.
    """
    print(f"\n[RIGID-FILE] Registering {file_source} → {file_target}")
    reg = pyregister.RigidFricpRegistration()
    reg.Paras_init(useinit=False, maxiter=100, stop=1e-5)
    reg.Reg(file_target, file_source, output_path)
    print(f"[RIGID-FILE] Completed → {output_path}")
    print("Deformed points shape:", reg.deformed_points_3X_.shape)


def run_rigid_numpy(target_pts, source_pts, target_n, source_n, output_path):
    """
    Perform rigid registration (FRICP method) using NumPy matrices.

    Parameters
    ----------
    target_pts : np.ndarray, shape (3, N)
        Target point cloud coordinates.
    source_pts : np.ndarray, shape (3, N)
        Source point cloud coordinates.
    target_n : np.ndarray, shape (3, N) or (0, 0)
        Target normals.
    source_n : np.ndarray, shape (3, N) or (0, 0)
        Source normals.
    output_path : str
        Output file path for the registered result.
    """
    print(f"\n[RIGID-NUMPY] Registering from NumPy matrices")
    reg = pyregister.RigidFricpRegistration()
    reg.Paras_init()
    reg.Read_data(target_pts, source_pts, target_n, source_n)
    reg.Register()
    reg.Output_data(output_path, "FRICP")
    print(f"[RIGID-NUMPY] Completed → {output_path}")
    print("Deformed points shape:", reg.deformed_points_3X_.shape)


def run_nonrigid_file(file_target, file_source, output_path):
    """
    Perform non-rigid registration (SPARE method) by directly reading files.

    Parameters
    ----------
    file_target : str
        Path to the target point cloud or mesh file.
    file_source : str
        Path to the source point cloud or mesh file.
    output_path : str
        Output file path for the registered result.
    """
    print(f"\n[NONRIGID-FILE] Registering {file_source} → {file_target}")
    reg = pyregister.NonrigidSpareRegistration()
    reg.Paras_init(iters = 30 ,stopcoarse=1e-3,stopfine=1e-4, uselandmark = False);
    reg.Reg(file_target, file_source, output_path)
    print(f"[NONRIGID-FILE] Completed → {output_path}")
    print("Deformed points shape:", reg.deformed_points_3X_.shape)



def run_nonrigid_numpy(target_pts, source_pts, target_n, source_n, output_path):
    """
    Perform non-rigid registration (SPARE method) using NumPy matrices.

    Parameters
    ----------
    target_pts : np.ndarray, shape (3, N)
        Target point cloud coordinates.
    source_pts : np.ndarray, shape (3, N)
        Source point cloud coordinates.
    target_n : np.ndarray, shape (3, N) or (0, 0)
        Target normals.
    source_n : np.ndarray, shape (3, N) or (0, 0)
        Source normals.
    output_path : str
        Output file path for the registered result.
    """
    print(f"\n[NONRIGID-NUMPY] Registering from NumPy matrices")
    reg = pyregister.NonrigidSpareRegistration()
    reg.Paras_init()
    reg.Read_data(target_pts, source_pts, target_n, source_n)
    reg.Register()
    reg.Output_data(output_path, "SPARE")
    print(f"[NONRIGID-NUMPY] Completed → {output_path}")
    print("Deformed points shape:", reg.deformed_points_3X_.shape)



def main():
    folder = "data"
    output_folder = os.path.join(folder, "results_register")
    os.makedirs(output_folder, exist_ok=True)
    target_ply = os.path.join(folder, "target.ply")

    source_ply = os.path.join(folder, "source.ply")

    target_obj = os.path.join(folder, "target.obj")

    source_obj = os.path.join(folder, "source.obj")

    # Rigid registration (file mode)
    if os.path.exists(target_ply) and os.path.exists(source_ply):
        run_rigid_file(target_ply, source_ply, os.path.join(output_folder, "rigid_file_result"))

        # Load into NumPy and run in-memory registration
        target_pts, target_n = read_vertices_normals(target_ply)
        source_pts, source_n = read_vertices_normals(source_ply)
        target_n = target_n if target_n is not None else np.empty((0, 0))
        source_n = source_n if source_n is not None else np.empty((0, 0))
        run_rigid_numpy(target_pts, source_pts, target_n, source_n,
                        os.path.join(output_folder, "rigid_numpy_result"))

    # Non-rigid registration (file mode)
    if os.path.exists( target_obj ) and os.path.exists(source_obj):
        run_nonrigid_file( target_obj , source_obj, os.path.join(output_folder, "nonrigid_file_result"))

        # Load into NumPy and run in-memory registration
        target_pts, target_n = read_vertices_normals( target_obj )
        source_pts, source_n = read_vertices_normals(source_obj)
        target_n = target_n if target_n is not None else np.empty((0, 0))
        source_n = source_n if source_n is not None else np.empty((0, 0))
        run_nonrigid_numpy(target_pts, source_pts, target_n, source_n,
                           os.path.join(output_folder, "nonrigid_numpy_result"))


if __name__ == "__main__":
    main()
