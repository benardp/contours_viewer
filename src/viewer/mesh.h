#ifndef MESH_H
#define MESH_H

#include <Eigen/Geometry>
#include <string>
#include <vector>

#include <surface_mesh.h>

#include "bvh.h"
#include "viewgraph.h"

class Subdiv;

class Mesh : public surface_mesh::Surface_mesh {

public:
  Mesh() : m_bvh(nullptr), m_nb_visible(-1) {}
  Mesh(const std::string &filename);

  ~Mesh();

  bool load(const std::string &filename);

  void init();

  const Vector3f get_mesh_center() { return m_bbox.center(); }

  real_t get_dist_max() { return m_bbox.diagonal().norm(); }

  /// returns face indices as a 3xN matrix of integers
  MapXu get_indices() { return MapXu(m_indices.data(), 3, n_faces()); }

  // Accessors to vertex attributes as Eigen's matrices:
  MapXf get_positions() {
    auto &vertices = get_vertex_property<Vector3f>("v:point").vector();
    return MapXf(vertices[0].data(), 3, vertices.size());
  }

  MapXf get_colors() {
    auto &vcolors = get_vertex_property<Vector3f>("v:color").vector();
    return MapXf(vcolors[0].data(), 3, vcolors.size());
  }

  MapXf get_normals() {
    auto &vnormals = get_vertex_property<Vector3f>("v:normal").vector();
    return MapXf(vnormals[0].data(), 3, vnormals.size());
  }

  /// Re-compute the aligned bounding box (needs to be called after editing
  /// vertex positions)
  void updateBoundingBox();

  const AABB &boundingBox() const { return m_bbox; }

  void tagConcaveEdges();

  void extractContours(ContourMode mode, const Vector3f &view_pos);

  void extractBoundaries();

  void extractBoundaryCurtainFolds(const Vector3f &view_pos);

  void extractSurfaceIntersections();

  void projectSVertices(const Matrix4f &model, const Matrix4f &proj,
                        const Vector2i &viewportSize);

  std::vector<Vector3f> &get_debug_points(nature_t nature);

  std::vector<Vector3f> &get_point_colors() { return m_point_colors; }

  const std::pair<int, int> &num2Dintersections() const {
    return m_2D_intersections;
  }

  void insertContours(const Vector3f &view_pos, const Subdiv *subdiv,
                      bool shift_allowed);

  void computeContourChains(ContourMode mode);

  void computeContourVisibility(const Vector3f &view_pos, ContourMode mode);

  std::vector<int> &get_chain_lengths() { return m_chain_lengths; }

  MapXu get_chain_indices() {
    return MapXu(m_chain_indices.data(), 1, m_chain_indices.size());
  }

  Matrix3Xf get_contours(VisibilityMode visible);

  Matrix3Xf get_chains(ContourMode mode);

  Matrix3Xf get_boundaries(bool visible_only);

  std::vector<int> &get_boundaries_lengths() { return m_boundaries_lengths; }

  Matrix3Xf &get_surface_intersections() { return m_surface_intersections; }

  std::vector<int> &get_surface_intersections_lengths() {
    return m_surface_intersections_lengths;
  }

  std::vector<std::list<Edge>> &get_chains() { return m_chains; }

  void computeSweepLineIntersections(ContourMode mode,
                                     const Matrix4f &projectionMatrix,
                                     const Matrix4f &viewMatrix,
                                     const Vector2i &viewport);

  bool hasBoundaries() const { return m_boundary_indices.size() > 0; }
  bool hasSurfaceIntersections() const {
    return m_surface_intersections.size() > 0;
  }

private:
  /// \returns true if d and e are on the same side of the abc triangle
  inline bool sameSide(const Eigen::Vector3d &a, const Eigen::Vector3d &b,
                       const Eigen::Vector3d &c, const Eigen::Vector3d &d,
                       const Eigen::Vector3d &e);

  // \returns returns a positive value if d is on the front side of triangle abc
  inline double frontSide(const Eigen::Vector3d &a, const Eigen::Vector3d &b,
                          const Eigen::Vector3d &c, const Eigen::Vector3d &d);

  /// \returns true if c and d are on the same side of the line ab
  inline bool sameSide(const Eigen::Vector2d &a, const Eigen::Vector2d &b,
                       const Eigen::Vector2d &c, const Eigen::Vector2d &d);

  void verticesOfFace(Face f, Vertex v[3]);
  void halfedgesOfFace(Face f, Halfedge h[3]);

  bool are_adjacent(Face f1, Face f2);

  AABB m_bbox;
  std::vector<uint32_t> m_indices;

  BVH *m_bvh;

  // Mesh contours
  std::vector<Edge> m_contour_edges;
  std::vector<std::list<Edge>> m_chains;
  std::vector<int> m_chain_lengths;
  std::vector<uint32_t> m_chain_indices;

  // Interpolated contours
  Matrix3Xf m_interpolated_contours;
  std::vector<int> m_interpolated_contour_faces;
  
  // Visibility
  std::vector<Vector3f> m_contour_sorted_by_visibility;
  int m_nb_visible;

  // Boundaries
  std::vector<uint32_t> m_boundary_indices;
  std::vector<int> m_boundaries_lengths;

  // Surface-Surface intersections
  Matrix3Xf m_surface_intersections;
  std::vector<int> m_surface_intersections_lengths;

  // View graph
  std::vector<SVertex *> m_svertices;
  std::vector<FEdge *> m_fedges;

  // Debug
  std::vector<Vector3f> ray_intersections;

  std::vector<Vector3f> m_debug_points, m_point_colors;
  std::pair<int, int> m_2D_intersections;
};

#endif // MESH_H
