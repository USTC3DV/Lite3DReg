<p align="right">
    <b> <img src="https://img.shields.io/badge/platforms-Windows%20%7C%20Linux-green" title="Supported Platforms"/> </b> <br>
<b>
  <img src="https://img.shields.io/badge/License-MIT-blue" title="License: MIT"/>
</b>
<br>

</p>

<img width="2064" height="384" alt="3" src="https://github.com/user-attachments/assets/cfa76e8e-cf73-44aa-a658-230e7e1ec5ef" />

#### Lightweight lib for 3D rigid registration and non-registration with robustness and efficiency
Lite3DRegLib, with Python and C++ APIs, is a lightweight 3D point cloud and mesh registration library supporting both rigid and non-rigid alignment. It is designed for research and educational purposes, while also serving as a solid foundation for developing advanced 3D applications. Compared to existing geometry and registration libraries (such as  [PCL ](https://pointclouds.org/),[Open3D](http://www.open3d.org/)), Lite3DRegLib is designed to be lightweight while using modern, high-performance, and robust registration algorithms, enabling fast and reliable 3D point cloud and mesh alignment.Designed to be accessible for beginners or users who only need registration results, Lite3DRegLib also includes interactive visualization and user-friendly tools via Gradio, making it easier to explore, analyze, and experiment with 3D data.
<a href="https://huggingface.co/spaces/USTC3DVer/Lite3DReg">
  <img src="https://img.shields.io/badge/%F0%9F%A4%97%20Gradio%20Demo-Huggingface-orange">
</a>


## Key features
- Lightweight:The source code is less than 0.3 MB.
  
- Efficient and Robust:By leveraging our algorithms, Lite3DRegLib significantly outperforms traditional registration methods such as ICP and CPD in terms of speed and robustness.Rigid registration algorithm:[Fast-Robust-ICP](https://github.com/yaoyx689/Fast-Robust-ICP);NonRigid registration algorithm:[spare](https://github.com/yaoyx689/spare).
  
- Easy and Flexible:Lite3DRegLib provides an interactive interface through Gradio, making it easy for users to visualize 3D registration results in real time. Users can conveniently adjust parameters and settings to explore different configurations, experiment with example datasets, and fine-tune the alignment process, all without writing additional code.




### Rigid registration [![IEEE 2021](https://img.shields.io/badge/IEEE_2021-00629B.svg?style=flat-square&logo=ieee&logoColor=white)](https://arxiv.org/abs/2007.07627)  

<table>
  <tr>
    <td width="50%">
      <img src="https://github.com/user-attachments/assets/ca99ec1f-3ca9-490f-8fe9-12ae599edc58" style="width:100%; height:auto;">
    </td>
    <td width="50%">
      <img src="https://github.com/user-attachments/assets/d0f332f5-80b9-4ea8-bc76-d19ca263e4d3" style="width:100%; height:auto;">
    </td>
  </tr>
</table>



### Non-rigid registration [![IEEE 2025](https://img.shields.io/badge/IEEE_2025-00629B.svg?style=flat-square&logo=ieee&logoColor=white)](https://arxiv.org/abs/2405.20188)

<table>
  <tr>
    <td>
      <img src="https://github.com/user-attachments/assets/220f4a79-221a-4350-86d9-2c633b56d6b2" width="400"/>
    </td>
    <td>
      <img src="https://github.com/user-attachments/assets/ffb851c7-a504-442b-8c08-419968203618" width="400"/>
    </td>
    <td>
      <img src="https://github.com/user-attachments/assets/7dd9fb53-0556-4879-9b5f-e6c8b2b83335" width="400"/>
    </td>
    <td>
      <img src="https://github.com/user-attachments/assets/fc1f11b6-8dab-45c1-ae70-8c1e6666ff34" width="400"/>
    </td>
  </tr>
</table>





## Huggingface & Gradio Visualization

Huggingface Space:
<a href="https://huggingface.co/spaces/USTC3DVer/Lite3DReg">
<img src="https://img.shields.io/badge/%F0%9F%A4%97%20Gradio%20Demo-Huggingface-orange">
</a>

Interactive Interface:

<img width="2094" height="1270" alt="image" src="https://github.com/user-attachments/assets/81d50639-ffd9-44a2-995a-3bebc664cc8c" />

You can adjust the registration parameters using the dropdown menus, and click the Reset to Default Parameters button to revert all settings to their default values.

## Lite3DReg repository layout
The repository includes a CMakeLists.txt file in its root directory, which serves as the main entry point for configuring and building programs, along with several subfolders.

- [3rd_party](./3rd_party) – source code of third-party libraries
- [cmake](./cmake) – CMake-related configuration files
- [cpp](./cpp) – source code of Lite3DReg
- [examples](./examples) – example data, C++ example, Python example
- [python](./python) – Python bindings for Lite3DReg
- [app.py](./app.py) – a Python visualization with Gradio
- [apt.txt](./apt.txt) – a list of required system packages
- [requirements.txt](./requirements.txt) – a list of required Python dependencies


## Dependencies

1. [Eigen-3.4.0](http://eigen.tuxfamily.org/index.php?title=Main_Page)
2. [OpenMesh-8.1](https://www.graphics.rwth-aachen.de/software/openmesh/)
3. (Optional for Linux) [MKL](https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl.html). It is used to speed up the LDLT solver(https://eigen.tuxfamily.org/dox-devel/group__PardisoSupport__Module.html). If it is not installed, Eigen's SimplicialLDLT will be used instead.

To install the required Python dependencies, use the following command:

    pip install numpy gradio plyfile trimesh plotly

## Compilation
Run `cd scripts && ./compile.sh` and a library and an executable  `Lightweight3DRegDemo.exe` of examples will be generated in the `build/examples` folder.

Or you can open CMakeLists.txt in Qtcreator, select the Release build mode, and compile.

## Interface
### C++
    RigidFricpRegistration *rigid;
    NonrigidSpareRegistration *spare;
    
    // Directly read data from files and output registration results
    rigid->Reg(target_file, source_file, outpath);
    spare->Reg(target_file, source_file, outpath);
    
    // Adjust parameters
    rigid->Paras_init(use_init=false, file_init="", max_iter=80, stop=1e-5);
    spare->Paras_init(int iters = 30 ,double stopcoarse=1e-3,double stopfine=1e-4,bool uselandmark = false,std::vector<int> src=std::vector<int>(),std::vector<int> tar=std::vector<int>());
    
    // Reset to default parameters
    rigid->Paras_init();
    spare->Paras_init();
    
    // Directly read vertices and normals (without using files)
    rigid->Read_data(target_points, source_points, target_normals, source_normals);
    spare->Read_data(target_points, source_points, target_normals, source_normals);
    
    // Perform registration
    rigid->Register();
    spare->Register();
    
    // Retrieve registration results
    Matrix3X result = rigid->deformed_points_3X_;
    Matrix3X result2 = spare->deformed_points_3X_;
    
    // Output results to files
    rigid->Output_data(outpath, "rigid");
    spare->Output_data(outpath, "spare");

### Python
    import pyregister
    
    # Create objects
    rigid = pyregister.RigidFricpRegistration()
    spare = pyregister.NonrigidSpareRegistration()
    
    # Directly read data from files and perform registration
    rigid.Reg("data/target.ply", "data/source.ply", "./outpath0/")
    spare.Reg("data/target.obj", "data/source.obj", "./outpath0/")
    
    # Optional parameter initialization
    rigid.Paras_init(useinit=False, fileinit="", maxiter=80, stop=1e-5)
    reg.Paras_init(iters = 30 ,stopcoarse=1e-3,stopfine=1e-4, uselandmark = False, src = [],tar = [])
    
    # Reset to default parameters
    rigid.Paras_init()
    spare.Paras_init()
    
    # Directly pass NumPy arrays of vertices and normals
    # target_p, source_p, target_n, source_n are all numpy.ndarray with shape=(3, n)
    rigid.Read_data(target_p, source_p, target_n, source_n)
    spare.Read_data(target_p, source_p, target_n, source_n)
    
    # Perform registration
    rigid.Register()
    spare.Register()
    
    # Retrieve registration results
    deformed_rigid = rigid.deformed_points_3X_  # numpy.ndarray, shape=(3, n)
    deformed_spare = spare.deformed_points_3X_
    
    # Output registration results to files
    rigid.Output_data("./outpath2/", "rigid")
    spare.Output_data("./outpath2/", "spare")

## parameters
### Rigid(fricp)
    Paras_init(use_init=false, file_init="", max_iter=80, stop=1e-5);
#### Initialization support 
If you have an initial transformation that can be applied on the input source model to roughly align with the input target model, you can set the first two parameters `.Paras_init(useinit=True
, fileinit="...")` to load your initial transformation. The format of the initial transformation is a 4x4 matrix(`[R, t; 0, 1]`), where `R` is a 3x3 rotation matrix and `t` is a 3x1 translation vector. These numbers are stored in 4 rows, and separated by spaces in each row. This format is the same as the output transformation of this code. It is worth mentioning that this code will align the center of gravity of the initial source and target models by default before starting the registration process, but this operation will be no longer used when the initial transformation is provided.
#### Iteration steps and Convergence accuracy
You can set the last two parameters of `Paras_init` to control the maximum number of iteration steps and the convergence tolerance. The registration process terminates either when the maximum number of iterations is reached or when the desired accuracy is achieved.
### NonRigid(spare)
    reg.Paras_init(iters = 30 ,stopcoarse=1e-3,stopfine=1e-4, uselandmark = False, src = [],tar = [])
#### Iteration steps and Convergence accuracy
Because our non-rigid registration consists of two stages, we define two convergence accuracy parameters — one for the coarse stage and another for the fine stage. Similar to the rigid registration process, each stage terminates either when the maximum number of iterations is reached or when the desired accuracy is achieved.
 
#### landmark point support

The registration process can optionally make use of landmark points to guide the deformation. When `uselandmark` is set to `True`, the algorithm incorporates corresponding landmark pairs from the source and target point sets (`src` and `tar`) into the registration. Here, `src` and `tar` represent the **IDs (indices)** of the corresponding points rather than their coordinates. These landmarks provide additional geometric constraints, improving the alignment stability and accuracy, especially in regions with sparse or ambiguous surface features.

## notes
The non-rigid registration method supports meshes or point clouds. When the input is a point cloud, the normal is required. When the source surface is represented as a point cloud, the deformation graph will be constructed by the farthest point sampling (FPS) based on Euclidean distance. 

## examples
The examples folder contains sample data and example programs. The steps to run the C++ example program are as follows:

1. Place main.cpp and the CMakeLists.txt from the C++ folder at the same level as the project folder.

1. In the CMakeLists.txt, replace 002hugging_face with the name of your project folder.

1. Compile using CMake.

Directory structure:

    workspace/
    │
    ├─ Lite3DRegLib/        # Library project folder
    │   ├─ cpp/             # Library source code (or rename to your project name)
    │   ├─ python/
    │   .
    │   .
    │   .
    │   └─ CMakeLists.txt
    │
    ├─ main.cpp             # Example program
    └─ CMakeLists.txt       # CMakeLists.txt for the example program; replace 002hugging_face with the project folder name (here it is Lite3DRegLib)


### Citation 
If you find our code or paper helps, please consider citing:

```
@article{zhang2021fast,
  author={Juyong Zhang and Yuxin Yao and Bailin Deng},
  title={Fast and Robust Iterative Closest Point}, 
  journal={IEEE Transactions on Pattern Analysis and Machine Intelligence}, 
  year={2022},
  volume={44},
  number={7},
  pages={3450-3466}}
```

```
@article{yao2025spare,
  author    = {Yao, Yuxin and Deng, Bailin and Hou, Junhui and Zhang, Juyong},
  title     = {SPARE: Symmetrized Point-to-Plane Distance for Robust Non-Rigid 3D Registration},
  journal   = {IEEE Transactions on Pattern Analysis and Machine Intelligence},
  year      = {2025},
  volume    = {},
  number    = {},
  pages     = {1-18},
}
```





