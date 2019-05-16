#ifndef COMMON_H
#define COMMON_H

#include <algorithm>
#include <atomic>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <queue>
#include <set>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>

#define GRAIN_SIZE 1024

// #define SINGLE_PRECISION
/* Application precision -- can be set to single or double precision */
#if defined(SINGLE_PRECISION)
typedef float real_t;
const float EPSILON = 1e-6;
#else
typedef double real_t;
const float EPSILON = 1e-8;
#endif

/* Useful Eigen typedefs based on the current precision */
typedef Eigen::Matrix<int32_t, 2, 1> Vector2i;
typedef Eigen::Matrix<int32_t, 3, 1> Vector3i;
typedef Eigen::Matrix<int32_t, 4, 1> Vector4i;
typedef Eigen::Matrix<uint32_t, 2, 1> Vector2u;
typedef Eigen::Matrix<uint32_t, 3, 1> Vector3u;
typedef Eigen::Matrix<uint32_t, 4, 1> Vector4u;
typedef Eigen::Matrix<uint8_t, 4, 1> Vector4u8;
typedef Eigen::Matrix<real_t, 2, 1> Vector2f;
typedef Eigen::Matrix<real_t, 3, 1> Vector3f;
typedef Eigen::Matrix<real_t, 4, 1> Vector4f;
typedef Eigen::Matrix<int32_t, Eigen::Dynamic, 1> VectorXi;
typedef Eigen::Matrix<uint32_t, Eigen::Dynamic, 1> VectorXu;
typedef Eigen::Matrix<uint8_t, Eigen::Dynamic, 1> VectorXu8;
typedef Eigen::Matrix<bool, Eigen::Dynamic, 1> VectorXb;
typedef Eigen::Matrix<real_t, Eigen::Dynamic, 1> VectorXf;
typedef Eigen::Matrix<int32_t, Eigen::Dynamic, Eigen::Dynamic> MatrixXi;
typedef Eigen::Matrix<uint32_t, Eigen::Dynamic, Eigen::Dynamic> MatrixXu;
typedef Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> MatrixXu8;
typedef Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> MatrixXf;
typedef Eigen::Matrix<real_t, 3, Eigen::Dynamic> Matrix3Xf;
typedef Eigen::Matrix<real_t, 2, Eigen::Dynamic> Matrix2Xf;
typedef Eigen::Matrix<real_t, 2, 2> Matrix2f;
typedef Eigen::Matrix<real_t, 3, 3> Matrix3f;
typedef Eigen::Matrix<real_t, 4, 4> Matrix4f;
typedef Eigen::Quaternion<real_t> Quaternionf;
typedef Eigen::Transform<real_t, 3, Eigen::Affine> Affine3f;
typedef Eigen::Translation<real_t, 3> Translation3f;
typedef Eigen::AngleAxis<real_t> AngleAxisf;

typedef Eigen::Map<MatrixXu> MapXu;
typedef Eigen::Map<MatrixXf> MapXf;

using std::cout;
using std::cerr;
using std::endl;

typedef enum { FRONT, BACK, CONTOUR, UNDEFINED } FacingType;
typedef enum { MESH_CONTOUR, INTERPOLATED_CONTOUR } ContourMode;
typedef enum { VISIBLE, INVISIBLE,  VISIBLE_AND_INVISIBLE } VisibilityMode;

typedef unsigned short nature_t;
namespace VertexNature {
static const nature_t S_VERTEX = 0;
static const nature_t INTERSECTION_3D = (1 << 1);
static const nature_t INTERSECTION_2D = (1 << 2);
static const nature_t BOUNDARY_CURTAIN_FOLD = (1 << 3);
static const nature_t CONTOUR_CURTAIN_FOLD = (1 << 4);
static const nature_t BIFURCATION = (1 << 5);
static const nature_t CURTAIN_FOLD =
    BOUNDARY_CURTAIN_FOLD | CONTOUR_CURTAIN_FOLD;
static const nature_t RAY_INTERSECTION = (1 << 6);
} // namespace VertexNature
namespace EdgeNature {
static const nature_t NO_FEATURE = 0;
static const nature_t SHARP_CONTOUR = (1 << 0);
static const nature_t SMOOTH_CONTOUR = (1 << 1);
static const nature_t CONTOUR = SHARP_CONTOUR | SMOOTH_CONTOUR;
static const nature_t BOUNDARY = (1 << 2);
static const nature_t SURFACE_INTERSECTION = (1 << 3);
static const nature_t IMAGE_INTERSECTION = (1 << 4);
} // namespace EdgeNature


typedef Vector4f Color;
const Color frontFaceColor = Color(254.0/255.0, 242.0/255.0, 192.0/255.0, 1.0);
const Color backFaceColor = Color(139.0/255.0, 191.0/255.0, 230.0/255.0, 1.0);
const Color undefinedFaceColor = Color(250.0/255.0, 50.0/255.0, 0.0, 1.0);
const Color inconsistentFaceColor = Color(0.8, 0.1, 0.8, 1.0);
const Color wireframeColor = Color(0.3, 0.3, 0.3, 1.0);
const Color normalsColor = Color(0.7, 0.7, 0., 1.0);
const Color contourColor = Color(0.6, 0.0, 0., 1.0);
const Color hiddenContourColor = Color(0.8, 0.8, 0.8, 1.0);
const Color boundariesColor = Color(0.1, 0.0, 0.8, 1.0);
const Color surfIntersectionsColor = Color(0.0, 0.5, 0.1, 1.0);
const Color intersection2DColor = Color(0.0, 1.0, 0.0, 1.0);
const Color intersection3DColor = Color(0.0, 1.0, 1.0, 1.0);
const Color curtainFoldColor = Color(1.0, 0.6, 0.0, 1.0);
const Color bifurcationColor = Color(1.0, 0.0, 0.0, 1.0);

// Paul Green-Armytage "A Colour Alphabet and the Limits of Colour Coding."
// (2010)
const Color alphabetColors[26] = {
    Color(240.0/255.0, 163.0/255.0, 255.0/255.0, 1.0), Color(0.0/255.0, 117.0/255.0, 220.0/255.0, 1.0),
    Color(153.0/255.0, 63.0/255.0, 0.0/255.0, 1.0),    Color(76.0/255.0, 0.0/255.0, 92.0/255.0, 1.0),
    Color(43.0/255.0, 206.0/255.0, 72.0/255.0, 1.0),   Color(255.0/255.0, 204.0/255.0, 153.0/255.0, 1.0),
    Color(128.0/255.0, 128.0/255.0, 128.0/255.0, 1.0), Color(148.0/255.0, 255.0/255.0, 181.0/255.0, 1.0),
    Color(143.0/255.0, 124.0/255.0, 0.0/255.0, 1.0),   Color(157.0/255.0, 204.0/255.0, 0.0/255.0, 1.0),
    Color(194.0/255.0, 0.0/255.0, 136.0/255.0, 1.0),   Color(0.0/255.0, 51.0/255.0, 128.0/255.0, 1.0),
    Color(25.0/255.0, 25.0/255.0, 25.0/255.0, 1.0),    Color(0.0/255.0, 92.0/255.0, 49.0/255.0, 1.0),
    Color(153.0/255.0, 0.0/255.0, 0.0/255.0, 1.0),     Color(255.0/255.0, 255.0/255.0, 128.0/255.0, 1.0),
    Color(255.0/255.0, 164.0/255.0, 5.0/255.0, 1.0),   Color(255.0/255.0, 168.0/255.0, 187.0/255.0, 1.0),
    Color(224.0/255.0, 255.0/255.0, 102.0/255.0, 1.0), Color(116.0/255.0, 10.0/255.0, 255.0/255.0, 1.0),
    Color(66.0/255.0, 102.0/255.0, 0.0/255.0, 1.0),    Color(255.0/255.0, 0.0/255.0, 16.0/255.0, 1.0),
    Color(94.0/255.0, 241.0/255.0, 242.0/255.0, 1.0),  Color(0.0/255.0, 153.0/255.0, 143.0/255.0, 1.0),
    Color(255.0/255.0, 255.0/255.0, 0.0/255.0, 1.0),   Color(255.0/255.0, 80.0/255.0, 5.0/255.0, 1.0)};

inline Vector3f project(const Vector3f &obj, 
                        const Matrix4f &model,
                        const Matrix4f &proj, 
                        const Vector2i &viewportSize) {
  Vector4f tmp;
  tmp << obj, 1;

  tmp = model * tmp;

  tmp = proj * tmp;

  tmp = tmp.array() / tmp(3);
  tmp = tmp.array() * 0.5f + 0.5f;
  tmp(0) = tmp(0) * viewportSize.x();
  tmp(1) = tmp(1) * viewportSize.y();

  return tmp.head(3);
}

inline Vector3f unproject(const Vector3f &win,
                          const Matrix4f &model,
                          const Matrix4f &proj,
                          const Vector2i &viewportSize) {
    Matrix4f Inverse = (proj * model).inverse();

    Vector4f tmp;
    tmp << win, 1;
    tmp(0) = tmp(0) / viewportSize.x();
    tmp(1) = tmp(1) / viewportSize.y();
    tmp = tmp.array() * 2.0f - 1.0f;

    Vector4f obj = Inverse * tmp;
    obj /= obj(3);

    return obj.head(3);
}

inline void projectToViewport(const Matrix3Xf &chain,
                              std::vector<Vector2f> &polyline,
                              const Matrix4f &projectionMatrix,
                              const Matrix4f &viewMatrix,
                              const Vector2i &viewport) {
  polyline.reserve(polyline.size() + chain.cols());
  for (int j = 0; j < chain.cols(); ++j) {
    Vector3f pos2D =
        project(chain.col(j), viewMatrix, projectionMatrix, viewport);
    polyline.push_back(pos2D.head<2>());
  }
}

template <typename T> void delete_pointed_to(T *const ptr) { delete ptr; }

inline std::string memString(size_t size, bool precise = false) {
  double value = (double)size;
  const char *suffixes[] = {"B", "KiB", "MiB", "GiB", "TiB", "PiB"};
  int suffix = 0;
  while (suffix < 5 && value > 1024.0f) {
    value /= 1024.0f;
    ++suffix;
  }

  std::ostringstream os;
  os << std::setprecision(suffix == 0 ? 0 : (precise ? 4 : 1)) << std::fixed
     << value << " " << suffixes[suffix];

  return os.str();
}

#endif