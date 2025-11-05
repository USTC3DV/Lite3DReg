#include "rigid_fricp_registration.h"
#include "nonrigid_spare_registration.h"  
#include "registration.h"
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>   // ✅ 支持 std::vector<int>
#include <Eigen/Dense>

namespace py = pybind11;

// ======================= 模块定义 ==========================
PYBIND11_MODULE(pyregister, m) {

    // -------- Registration 基类 --------
    py::class_<Registration>(m, "Registration")
        .def("Read_data",
             py::overload_cast<const std::string&, const std::string&>(&Registration::Read_data),
             py::arg("file_target"),
             py::arg("file_source"))
        .def("Read_data",
             [](Registration &self,
                const Eigen::Ref<const Eigen::MatrixXd> &target_p,
                const Eigen::Ref<const Eigen::MatrixXd> &source_p,
                const Eigen::Ref<const Eigen::MatrixXd> &target_n,
                const Eigen::Ref<const Eigen::MatrixXd> &source_n)
             {
                 std::cout << "[pybind] entering Registration::Read_data" << std::endl;

                 Matrix3X tp = target_p.topRows(3);
                 Matrix3X sp = source_p.topRows(3);
                 Matrix3X tn = target_n.topRows(3);
                 Matrix3X sn = source_n.topRows(3);

                 self.Read_data(tp, sp, tn, sn);

                 std::cout << "[pybind] leaving Registration::Read_data" << std::endl;
             },
             py::arg("target_p"),
             py::arg("source_p"),
             py::arg("target_n") = Eigen::MatrixXd(),
             py::arg("source_n") = Eigen::MatrixXd()
             )


        .def("Register", &Registration::Register)
        .def("Output_data", &Registration::Output_data,
             py::arg("out_path"),
             py::arg("method_name"));

    // -------- 暴露 spare::Parameters --------
    py::class_<spare::Parameters>(m, "SpareParameters")
        .def(py::init<>())
        .def_readwrite("use_landmark", &spare::Parameters::use_landmark)
        .def_readwrite("landmark_src", &spare::Parameters::landmark_src)
        .def_readwrite("landmark_tar", &spare::Parameters::landmark_tar);

    // -------- Rigid FRICP Registration --------
    py::class_<RigidFricpRegistration, Registration>(m, "RigidFricpRegistration")
        .def(py::init<>())
        .def("Reg", &RigidFricpRegistration::Reg,
             py::arg("file_target"),
             py::arg("file_source"),
             py::arg("out_path"))
        .def("Paras_init", &RigidFricpRegistration::Paras_init,
             py::arg("useinit") = false,
             py::arg("fileinit") = " ",
             py::arg("maxiter") = 100,
             py::arg("stop") = 1e-5)
        .def_readwrite("deformed_points_3X_", &RigidFricpRegistration::deformed_points_3X_);

    // -------- Nonrigid Spare Registration --------
    py::class_<NonrigidSpareRegistration, Registration>(m, "NonrigidSpareRegistration")
        .def(py::init<>())
        .def("Reg", &NonrigidSpareRegistration::Reg,
             py::arg("file_target"),
             py::arg("file_source"),
             py::arg("out_path"))
        .def("Paras_init", &NonrigidSpareRegistration::Paras_init,
     py::arg("iters") = 30,
     py::arg("stopcoarse") = 1e-3,
     py::arg("stopfine") = 1e-4,
     py::arg("uselandmark") = false,
     py::arg("src") = std::vector<int>(),
     py::arg("tar") = std::vector<int>())

        .def_readwrite("deformed_points_3X_", &NonrigidSpareRegistration::deformed_points_3X_)
        .def_readwrite("paras", &NonrigidSpareRegistration::paras);  // ✅ 暴露 paras
}
