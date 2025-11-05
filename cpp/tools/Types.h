#ifndef TYPES_H
#define TYPES_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>

#include <iostream>
#include <fstream>
#include <iomanip>

#include <stdlib.h>
#include <stdio.h>
#include <cstdio>

#include <algorithm>
#include <math.h>
#include <cmath>

#include <vector>
#include <list>
#include <set>
#include <queue>
#include <map>
#include <string>

#include <time.h>
#include <assert.h>
#include <cstddef>
#include <limits>

//#define USE_FLOAT_SCALAR

#define SAME_THRESHOLD 1e-6


#ifdef USE_FLOAT_SCALAR
typedef float Scalar;
#else
typedef double Scalar;
#endif

typedef Eigen::SparseMatrix<Scalar, Eigen::ColMajor> ColMajorSparseMatrix;
typedef Eigen::SparseMatrix<Scalar, Eigen::RowMajor> RowMajorSparseMatrix;


#ifdef EIGEN_DONT_ALIGN
#define EIGEN_ALIGNMENT Eigen::DontAlign
#else
#define EIGEN_ALIGNMENT Eigen::AutoAlign
#endif


typedef Eigen::Triplet<Scalar> Triplet;

// ======================================================
// 基础模板别名
// ======================================================
template <int Rows, int Cols, int Options = (Eigen::ColMajor | Eigen::AutoAlign)>
using MatrixT = Eigen::Matrix<Scalar, Rows, Cols, Options>;

// ======================
// 向量类型
// ======================
using Vector2 = MatrixT<2,1>;
using Vector3 = MatrixT<3,1>;
using Vector4 = MatrixT<4,1>;
using Vector6 = MatrixT<6,1>;
using VectorX = MatrixT<Eigen::Dynamic,1>;
using VectorN = Vector3;                // 别名（等价于 Vector3）

// ======================
// 固定尺寸矩阵
// ======================
using Matrix22 = MatrixT<2,2>;
using Matrix23 = MatrixT<2,3>;
using Matrix32 = MatrixT<3,2>;
using Matrix33 = MatrixT<3,3>;
using Matrix34 = MatrixT<3,4>;
using Matrix44 = MatrixT<4,4>;
using Matrix66 = MatrixT<6,6>;
using EigenMatrix12 = MatrixT<12,12>;

// ======================
// 动态尺寸矩阵
// ======================
using Matrix2X = MatrixT<2,Eigen::Dynamic>;
using Matrix3X = MatrixT<3,Eigen::Dynamic>;
using Matrix4X = MatrixT<4,Eigen::Dynamic>;
using MatrixX2 = MatrixT<Eigen::Dynamic,2>;
using MatrixX3 = MatrixT<Eigen::Dynamic,3>;
using MatrixXX = MatrixT<Eigen::Dynamic,Eigen::Dynamic>;
using Vertices = Matrix3X;              // 顶点矩阵，3×N

// ======================
// 特殊变换与几何类型
// ======================
using Affine3  = Eigen::Transform<Scalar,3,Eigen::Affine>;
using Affine3d = Affine3;               // 同义，方便兼容
using AffineMatrix3 = Matrix44;         // 4×4 变换矩阵

using EigenAngleAxis  = Eigen::AngleAxis<Scalar>;
using EigenQuaternion = Eigen::Quaternion<Scalar, Eigen::DontAlign>;

// ======================
// 块类型
// ======================
using Block33 = Eigen::Block<Matrix66,3,3>;

// Conversion between a 3d vector type to Eigen::Vector3d
template<typename Vec_T>
inline Vector3 to_eigen_vec3(const Vec_T &vec)
{
    return Vector3(vec[0], vec[1], vec[2]);
}


template<typename Vec_T>
inline Vec_T from_eigen_vec3(const Vector3 &vec)
{
    Vec_T v;
    v[0] = vec(0);
    v[1] = vec(1);
    v[2] = vec(2);

    return v;
}

enum VertexState
{
    OUTSIDE,
    FRONT,
    INSIDE
};

struct TriTraits : public OpenMesh::DefaultTraits {
    #ifdef USE_FLOAT_SCALAR
    typedef OpenMesh::Vec3f Point;
    typedef OpenMesh::Vec3f Normal;
#else
    typedef OpenMesh::Vec3d Point;
    typedef OpenMesh::Vec3d Normal;
#endif
    typedef OpenMesh::Vec4f Color;



    VertexAttributes(OpenMesh::Attributes::Status);
    FaceAttributes(OpenMesh::Attributes::Status);
    EdgeAttributes(OpenMesh::Attributes::Status);
    HalfedgeAttributes(OpenMesh::Attributes::Status);

    VertexTraits
    {
        VertexT() : geodesic_distance(1e100), state(OUTSIDE),incident_point(0.0),
          saddle_or_boundary(false)
        {};
    public:
        Scalar geodesic_distance;
        VertexState state;
        OpenMesh::FaceHandle incident_face;
        Scalar incident_point;
        bool saddle_or_boundary;
    };
    FaceTraits
    {
        FaceT() : corner_angles(0,0,0)
        {};
        public:
        Vector3 corner_angles;
    };
    EdgeTraits
    {
        EdgeT() : length(0.0)
        {};
    public:
        Scalar length;
    };
};
typedef OpenMesh::TriMesh_ArrayKernelT<TriTraits> Mesh;
#ifdef USE_FLOAT_SCALAR
typedef OpenMesh::Vec3f Vec3;
#else
typedef OpenMesh::Vec3d Vec3;
#endif


class Matrix3333 // 3x3 matrix: each element is a 3x3 matrix
{
public:
    Matrix3333();
    Matrix3333(const Matrix3333& other);
    ~Matrix3333() {}

    void SetZero(); // [0 0 0; 0 0 0; 0 0 0]; 0 = 3x3 zeros
    void SetIdentity(); //[I 0 0; 0 I 0; 0 0 I]; 0 = 3x3 zeros, I = 3x3 identity

                        // operators
    Matrix33& operator() (int row, int col);
    Matrix3333 operator+ (const Matrix3333& plus);
    Matrix3333 operator- (const Matrix3333& minus);
    Matrix3333 operator* (const Matrix33& multi);
    friend Matrix3333 operator* (const Matrix33& multi1, Matrix3333& multi2);
    Matrix3333 operator* (Scalar multi);
    friend Matrix3333 operator* (Scalar multi1, Matrix3333& multi2);
    Matrix3333 transpose();
    Matrix33 Contract(const Matrix33& multi); // this operator is commutative
    Matrix3333 Contract(Matrix3333& multi);

    //protected:

    Matrix33 mat[3][3];
};

class Matrix2222 // 2x2 matrix: each element is a 2x2 matrix
{
public:
    Matrix2222();
    Matrix2222(const Matrix2222& other);
    ~Matrix2222() {}

    void SetZero(); // [0 0; 0 0]; 0 = 2x2 zeros
    void SetIdentity(); //[I 0; 0 I;]; 0 = 2x2 zeros, I = 2x2 identity

                        // operators and basic functions
    Matrix22& operator() (int row, int col);
    Matrix2222 operator+ (const Matrix2222& plus);
    Matrix2222 operator- (const Matrix2222& minus);
    Matrix2222 operator* (const Matrix22& multi);
    friend Matrix2222 operator* (const Matrix22& multi1, Matrix2222& multi2);
    Matrix2222 operator* (Scalar multi);
    friend Matrix2222 operator* (Scalar multi1, Matrix2222& multi2);
    Matrix2222 transpose();
    Matrix22 Contract(const Matrix22& multi); // this operator is commutative
    Matrix2222 Contract(Matrix2222& multi);

protected:

    Matrix22 mat[2][2];
};

// dst = src1 \kron src2
void directProduct(Matrix3333& dst, const Matrix33& src1, const Matrix33& src2);
void directProduct(Matrix2222& dst, const Matrix22& src1, const Matrix22& src2);
#endif // TYPES_H
///////////////////////////////////////////////////////////////////////////////
