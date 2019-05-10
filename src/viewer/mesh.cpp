#include "mesh.h"
#include "svg.h"
#include "sweepLine.h"

#include <geogram/numerics/predicates.h>

Mesh::Mesh(const std::string &filename) : m_bvh(nullptr), m_visible(-1) {
  load(filename);
}

Mesh::~Mesh() {
  std::for_each(m_fedges.begin(), m_fedges.end(), delete_pointed_to<FEdge>);
  m_fedges.clear();
  std::for_each(m_svertices.begin(), m_svertices.end(),
                delete_pointed_to<SVertex>);
  m_svertices.clear();
}

void Mesh::init() {
  if (!is_triangle_mesh()) {
    triangulate();
  }

  updateBoundingBox();

  update_face_normals();

  auto vpositions = get_vertex_property<Vector3f>("v:point");
  assert(vpositions);

  auto vnormals = get_vertex_property<Vector3f>("v:normal");
  if (!vnormals) {
    update_vertex_normals();
    vnormals = get_vertex_property<Vector3f>("v:normal");
  }
  assert(vnormals);

  auto colors = get_vertex_property<Vector3f>("v:color");
  if (!colors) {
    colors = add_vertex_property<Vector3f>("v:color");
    this->get_colors().setOnes();
  }

  m_indices.reserve(3 * n_faces());
  // face iterator
  Face_iterator fit, fend = faces_end();
  // vertex circulator
  Vertex vertices[3];
  for (fit = faces_begin(); fit != fend; ++fit) {
    verticesOfFace(*fit, vertices);
    for (Vertex v : vertices)
      m_indices.push_back(v.idx());
  }
}

bool Mesh::load(const std::string &filename) {

  if (!read(filename))
    return false;

  init();

  return true;
}

void Mesh::updateBoundingBox() {
  m_bbox.clear();
  auto vertices = get_vertex_property<Vector3f>("v:point");
  for (const auto &p : vertices.vector())
    m_bbox.expandBy(p);
}

void Mesh::tagConcaveEdges() {
  auto concave = edge_property<bool>("e:concave");
  concave.vector().assign(concave.vector().size(), false);
  auto vpositions = get_vertex_property<Vector3f>("v:point");

  for (uint32_t i = 0; i != n_edges(); ++i) {
    Edge e = Edge(i);
    Halfedge h = halfedge(e, 0);

    const Vector3f &a = vpositions[from_vertex(h)];
    const Vector3f &b = vpositions[to_vertex(h)];

    // concavity test
    const Vector3f &c = vpositions[to_vertex(next_halfedge(h))];
    const Vector3f &d =
        vpositions[to_vertex(next_halfedge(opposite_halfedge(h)))];
    if (frontSide(a.cast<double>(), b.cast<double>(), c.cast<double>(),
                  d.cast<double>()) > 0) {
      // concave edge
      concave[e] = true;
    }
  }
}

inline real_t find_zero_linear(real_t val0, real_t val1) {
  return val0 / (val0 - val1);
}

inline void interpolate(real_t d0, real_t d1, real_t d2, const Vector3f &v0,
                        const Vector3f &v1, const Vector3f &v2, Vector3f &vm0,
                        Vector3f &vm1) {
  real_t w0 = find_zero_linear(d0, d1);
  real_t w1 = find_zero_linear(d0, d2);
  vm0 = w0 * v1 + (1.f - w0) * v0;
  vm1 = w1 * v2 + (1.f - w1) * v0;
  if (d0 > 0) {
    std::swap(vm0, vm1);
  }
}

template <typename T> inline int sign(T val) {
  return (T(0) < val) - (val < T(0));
}

void Mesh::extractContours(ContourMode mode, const Vector3f &view_pos) {
  auto vpositions = get_vertex_property<Vector3f>("v:point");
  auto is_contour = get_edge_property<real_t>("e:contour");
  if (!is_contour) {
    is_contour = add_edge_property<real_t>("e:contour", -1.f);
  } else {
    for (size_t i = 0; i < is_contour.vector().size(); ++i)
      is_contour.vector()[i] = -1.f;
  }

  m_contour_edges.clear();
  m_contour_edges.reserve(std::sqrt(n_faces()));

  auto fnormals = get_face_property<Vector3f>("f:normal");
  auto ndotv = vertex_property<real_t>("v:ndotv");

  // compute n dot v
  if (mode == INTERPOLATED_CONTOUR) {
    auto vnormals = get_vertex_property<Vector3f>("v:normal");
    for (uint32_t i = 0; i != n_vertices(); ++i) {
      Vertex v = Vertex(i);
      ndotv[v] = vnormals[v].dot((view_pos - vpositions[v]).normalized());
    }
  }

  for (uint32_t i = 0; i < n_edges(); ++i) {
    Edge e = Edge(i);

    if (mode == INTERPOLATED_CONTOUR) {
      Vertex i0 = vertex(e, 0);
      Vertex i1 = vertex(e, 1);
      real_t d0 = ndotv[i0];
      real_t d1 = ndotv[i1];
      if ((d0 > 0 && d1 <= 0) || (d0 < 0 && d1 >= 0)) {
        real_t w = find_zero_linear(d0, d1);
        m_contour_edges.push_back(e);
        is_contour[e] = w;
      }
    } else {
      if (is_boundary(e))
        continue;
      Face f1 = face(e, 0);
      Face f2 = face(e, 1);
      const Vector3f &normal1 = fnormals[f1];
      const Vector3f &normal2 = fnormals[f2];
      const Vector3f &v1 = vpositions[vertex(e, 0)];
      const Vector3f &v2 = vpositions[vertex(e, 1)];

      if (sign((view_pos - v1).dot(normal1)) !=
          sign((view_pos - v2).dot(normal2))) {
        // silhouette edge found
        m_contour_edges.push_back(e);
        is_contour[e] = 1.f;
      }
    }
  }

  if (mode == MESH_CONTOUR) {
    m_mesh_contours_indices = MatrixXu(2, m_contour_edges.size());
    for (size_t i = 0; i < m_contour_edges.size(); ++i) {
      m_mesh_contours_indices.col(i) =
          Vector2u(vertex(m_contour_edges.at(i), 0).idx(),
                   vertex(m_contour_edges.at(i), 1).idx());
    }
  } else {
    std::vector<std::tuple<Vector3f, Vector3f, int>> contour_pairs;

    Halfedge_around_face_circulator hit, hit_end;
    for (uint32_t i = 0; i < n_faces(); ++i) {
      Face f = Face(i);
      hit = hit_end = halfedges(f);
      Edge e[2];
      real_t w[2];
      int j = 0;
      do {
        if (is_contour[edge(*hit)] >= 0) {
          e[j] = edge(*hit);
          w[j] = is_contour[edge(*hit)];
          j++;
        }
      } while (++hit != hit_end && j < 2);
      assert(e[0].is_valid() == e[1].is_valid());
      if (!e[0].is_valid() || !e[1].is_valid())
        continue;
      Vector3f vm[2];
      for (j = 0; j < 2; j++) {
        Vertex i0 = vertex(e[j], 0);
        Vertex i1 = vertex(e[j], 1);
        const Vector3f &v0 = vpositions[i0];
        const Vector3f &v1 = vpositions[i1];
        vm[j] = w[j] * v1 + (1.f - w[j]) * v0;
      }
      contour_pairs.push_back(std::make_tuple(vm[0], vm[1], i));
    }

    m_interpolated_contours = Matrix3Xf(3, contour_pairs.size() * 2);
    m_interpolated_contour_faces.resize(contour_pairs.size());
    for (size_t i = 0; i < contour_pairs.size(); ++i) {
      m_interpolated_contours.col(i * 2) = std::get<0>(contour_pairs.at(i));
      m_interpolated_contours.col(i * 2 + 1) = std::get<1>(contour_pairs.at(i));
      m_interpolated_contour_faces[i] = std::get<2>(contour_pairs.at(i));
    }
  }
}

void Mesh::extractBoundaries() {
  m_boundary_indices.clear();
  auto visited = get_edge_property<real_t>("e:visited");
  if (!visited) {
    visited = add_edge_property<real_t>("e:visited", false);
  } else {
    for (size_t i = 0; i < visited.vector().size(); ++i)
      visited.vector()[i] = false;
  }
  auto boundary_fedges = edge_property<FEdge *>("e:fedge", nullptr);
  auto svertices = vertex_property<SVertex *>("v:svertex", nullptr);
  auto vpositions = get_vertex_property<Vector3f>("v:point");

  Halfedge_around_vertex_circulator hit, hit_end;
  int size = 0;
  for (Edge e : edges()) {
    if (is_boundary(e) && !visited[e]) {
      Halfedge sh = halfedge(e, 0);
      if (is_boundary(sh))
        sh = halfedge(e, 1);
      m_boundary_indices.push_back(from_vertex(sh).idx());
      const Vector3f &p = vpositions[from_vertex(sh)];
      SVertex *sv_first, *sv0;
      sv_first = new SVertex(p);
      sv0 = sv_first;
      assert(svertices[from_vertex(sh)] == nullptr);
      svertices[from_vertex(sh)] = sv0;
      m_svertices.push_back(sv0);
      visited[e] = true;
      Halfedge h = sh;
      do {
        hit = hit_end = halfedges(to_vertex(h));
        do {
          if (edge(h) != edge(*hit) && is_boundary(edge(*hit))) {
            h = *hit;
            break;
          }
        } while (++hit != hit_end);
        assert(hit != hit_end);
        visited[edge(h)] = true;
        SVertex *sv1;
        if (svertices[from_vertex(h)]) {
          sv1 = svertices[from_vertex(h)];
        } else {
          const Vector3f &p = vpositions[from_vertex(h)];
          sv1 = new SVertex(p);
          svertices[from_vertex(h)] = sv1;
          m_svertices.push_back(sv1);
        }
        m_fedges.push_back(
            new FEdgeSharp(sv0, sv1, edge(h), EdgeNature::BOUNDARY));
        boundary_fedges[edge(h)] = m_fedges.back();
        sv0->m_fedges.push_back(m_fedges.back());
        sv1->m_fedges.push_back(m_fedges.back());
        sv0 = sv1;
        m_boundary_indices.push_back(from_vertex(h).idx());
      } while (edge(h) != edge(sh));
      m_boundaries_lengths.push_back(m_boundary_indices.size() - size);
      size = m_boundary_indices.size();
    }
  }
}

void Mesh::extractBoundaryCurtainFolds(const Vector3f &c) {
  auto vpositions = get_vertex_property<Vector3f>("v:point");
  auto svertices = get_vertex_property<SVertex *>("v:svertex");
  assert(svertices);

  for (uint32_t i = 0; i != m_fedges.size(); ++i) {
    FEdge *fe = m_fedges[i];
    if (!fe || fe->m_nature != EdgeNature::BOUNDARY)
      continue;
    FEdgeSharp *fs = dynamic_cast<FEdgeSharp *>(fe);
    assert(fs);
    Edge edge = fs->m_edge;
    for (int i = 0; i < 2; ++i) {
      Vertex ie = vertex(edge, i);
      const Vector3f &e = vpositions[ie];
      Vertex ip = vertex(edge, (i + 1) % 2);
      const Vector3f &p = vpositions[ip];
      Halfedge_around_vertex_circulator hit, hit_end;
      hit = hit_end = halfedges(ip);
      do {
        if (is_boundary(*hit))
          continue;
        Vertex iq = to_vertex(*hit);
        Vertex ir = to_vertex(next_halfedge(*hit));
        if (iq == ie || ir == ie || iq == ip || ir == ip)
          continue;
        const Vector3f &q = vpositions[iq];
        const Vector3f &r = vpositions[ir];
        if (!sameSide(p.cast<double>(), q.cast<double>(), r.cast<double>(),
                      c.cast<double>(), e.cast<double>()) &&
            sameSide(c.cast<double>(), p.cast<double>(), q.cast<double>(),
                     e.cast<double>(), r.cast<double>()) &&
            sameSide(c.cast<double>(), p.cast<double>(), r.cast<double>(),
                     e.cast<double>(), q.cast<double>())) {
          // boundary curtain fold found
          SVertex *sv = svertices[ip];
          assert(sv);
          sv->m_nature = VertexNature::BOUNDARY_CURTAIN_FOLD | sv->m_nature;
        }
      } while (++hit != hit_end);
    }
  }
}

extern int tri_tri_intersection_test_3d(double p1[3], double q1[3],
                                        double r1[3], double p2[3],
                                        double q2[3], double r2[3],
                                        int *coplanar, double source[3],
                                        double target[3]);

void Mesh::verticesOfFace(Face f, Vertex v[3]) {
  Vertex_around_face_circulator fvit, fvend;
  fvit = fvend = vertices(f);
  v[0] = *fvit;
  ++fvit;
  v[2] = *fvit;
  do {
    v[1] = v[2];
    ++fvit;
    v[2] = *fvit;
  } while (++fvit != fvend);
}

void Mesh::halfedgesOfFace(Face f, Halfedge h[3]) {
  Halfedge_around_face_circulator fvit, fvend;
  fvit = fvend = halfedges(f);
  h[0] = *fvit;
  ++fvit;
  h[2] = *fvit;
  do {
    h[1] = h[2];
    ++fvit;
    h[2] = *fvit;
  } while (++fvit != fvend);
}

bool Mesh::are_adjacent(Face f1, Face f2) {
  Vertex_around_face_circulator vit, vit_end;
  vit = vit_end = vertices(f2);
  do {
    Face_around_vertex_circulator fit, fit_end;
    fit = fit_end = faces(*vit);
    do {
      if (*fit == f1)
        return true;
    } while (++fit != fit_end);
  } while (++vit != vit_end);
  return false;
};

void Mesh::extractSurfaceIntersections() {

  if (m_bvh)
    delete m_bvh;
  m_bvh = new BVH(
      &m_indices, &get_vertex_property<Vector3f>("v:point").vector(),
      &get_vertex_property<Vector3f>("v:normal").vector(), boundingBox());
  m_bvh->build();

  auto vpositions = get_vertex_property<Vector3f>("v:point");
  struct SurfaceIntersection {
    SurfaceIntersection(real_t _t, Edge _e, const Vector3f &_p)
        : t(_t), e(_e), visited(false) {
      neighboors.reserve(2);
      sv = new SVertex(_p);
    }
    real_t t;
    Edge e;
    bool visited;
    std::vector<SurfaceIntersection *> neighboors;
    SVertex *sv;
  };
  auto surf_intersect =
      add_edge_property<std::vector<SurfaceIntersection *>>("e:intersection");
  auto boundary_fedges = edge_property<FEdge *>("e:fedge", nullptr);
  auto fedges_f = face_property<std::vector<FEdgeIntersection *>>("f:fedge");
  auto svertices = edge_property<std::vector<SVertex *>>("e:svertex");

  std::set<Edge> intersected_edges;

  auto findSurfaceIntersection = [&](Edge e,
                                     real_t t) -> SurfaceIntersection * {
    SurfaceIntersection *s = nullptr;
    for (size_t j = 0; j < surf_intersect[e].size(); j++) {
      if (std::abs(surf_intersect[e][j]->t - t) < EPSILON) {
        // same intersection point
        s = surf_intersect[e][j];
        break;
      }
    }
    if (!s) {
      const Vector3f &v0 = vpositions[vertex(e, 0)];
      const Vector3f &v1 = vpositions[vertex(e, 1)];
      Vector3f p = v0 + t * (v1 - v0).normalized();
      s = new SurfaceIntersection(t, e, p);
      surf_intersect[e].push_back(s);
      svertices[e].push_back(s->sv);
      m_svertices.push_back(s->sv);

      if (boundary_fedges[e]) { // 3D intersection with a boundary edge
        s->sv->m_nature = VertexNature::INTERSECTION_3D | s->sv->m_nature;
      }
    }
    return s;
  };

  for (Face f1 : faces()) {
    Vertex v1[3];
    verticesOfFace(f1, v1);
    Eigen::Vector3d p1 = vpositions[v1[0]].cast<double>();
    Eigen::Vector3d q1 = vpositions[v1[1]].cast<double>();
    Eigen::Vector3d r1 = vpositions[v1[2]].cast<double>();

    auto intersect = [&](const Eigen::Vector3d &o, const Eigen::Vector3d &d) {
      Ray r = Ray(o.cast<real_t>(), d.normalized().cast<real_t>(), 0, d.norm());
      std::vector<uint32_t> indices;
      if (m_bvh->rayIntersect(r, indices)) {
        for (int idx : indices) {
          Face f2 = Face(idx);

          if (f1 == f2 || are_adjacent(f1, f2))
            continue;
          bool found = false;
          for (FEdgeIntersection *fe : fedges_f[f1]) {
            if (fe->m_f1 == f2 || fe->m_f2 == f2) {
              found = true;
              break;
            }
          }
          if (found) // Face already processed
            continue;

          Vertex v2[3];
          verticesOfFace(f2, v2);
          Eigen::Vector3d p2 = vpositions[v2[0]].cast<double>();
          Eigen::Vector3d q2 = vpositions[v2[1]].cast<double>();
          Eigen::Vector3d r2 = vpositions[v2[2]].cast<double>();
          Eigen::Vector3d i1, i2;
          int colplanar;
          int res = tri_tri_intersection_test_3d(
              p1.data(), q1.data(), r1.data(), p2.data(), q2.data(), r2.data(),
              &colplanar, i1.data(), i2.data());

          if (res == 0) // no intersection
            continue;

          // figure out which intersection goes with which face edge
          Face faces[2] = {f1, f2};
          real_t min_distA, min_distB;
          min_distA = min_distB = std::numeric_limits<real_t>::infinity();
          Edge eA, eB;
          real_t tA, tB;
          for (int i = 0; i < 2; i++) {
            Halfedge_around_face_circulator hit, hit_end;
            hit = hit_end = halfedges(faces[i]);
            do {
              Edge e = edge(*hit);
              const Vector3f &v0 = vpositions[vertex(e, 0)];
              const Vector3f &v1 = vpositions[vertex(e, 1)];
              Vector3f dir = (v1 - v0).normalized();
              Vector3f e1 = v0 - i1.cast<real_t>();
              real_t dist = dir.cross(e1).squaredNorm();
              if (dist <= min_distA) {
                min_distA = dist;
                tA = -dir.dot(e1);
                eA = e;
              }
              Vector3f e2 = v0 - i2.cast<real_t>();
              dist = dir.cross(e2).squaredNorm();
              if (dist <= min_distB) {
                min_distB = dist;
                tB = -dir.dot(e2);
                eB = e;
              }
            } while (++hit != hit_end);
          }

          if (!eA.is_valid() || !eB.is_valid())
            throw std::runtime_error(
                "Cannot find matching surface-surface intersections.");

          SurfaceIntersection *sA = findSurfaceIntersection(eA, tA);
          SurfaceIntersection *sB = findSurfaceIntersection(eB, tB);
          sA->neighboors.push_back(sB);
          sB->neighboors.push_back(sA);

          intersected_edges.insert(eA);
          intersected_edges.insert(eB);

          FEdgeIntersection *fe = new FEdgeIntersection(sA->sv, sB->sv, f1, f2);
          m_fedges.push_back(fe);
          sA->sv->m_fedges.push_back(fe);
          sB->sv->m_fedges.push_back(fe);
          fedges_f[f1].push_back(fe);
          fedges_f[f2].push_back(fe);
        }
      }
    };

    intersect(p1, q1 - p1);
    intersect(p1, r1 - p1);
    intersect(q1, r1 - q1);
  }

  // Chaining
  std::vector<std::list<SVertex *>> intersections_chains;
  int sum = 0;
  for (auto it = intersected_edges.begin(); it != intersected_edges.end();
       it++) {
    Edge e = *it;
    for (SurfaceIntersection *s : surf_intersect[e]) {
      if (s->visited)
        continue;
      SurfaceIntersection *is = s;
      std::list<SVertex *> chain;
      bool forward = true;
      while (s) {
        s->visited = true;

        if (forward)
          chain.push_back(s->sv);
        else
          chain.push_front(s->sv);

        bool found = false;
        bool loop_found = false;
        for (auto ns : s->neighboors) {
          if (!ns->visited) {
            s = ns;
            loop_found = false;
            found = true;
            break;
          }
          if (ns == is)
            loop_found = true;
        }
        if (loop_found) {
          chain.push_front(s->sv);
          break;
        }
        if (!found) {
          if (forward) {
            for (auto ns : is->neighboors) {
              if (!ns->visited) {
                s = ns;
                forward = false;
                found = true;
                break;
              }
            }
          }
          if (!found)
            s = nullptr;
        }
      }
      intersections_chains.push_back(chain);
      m_surface_intersections_lengths.push_back(chain.size());
      sum += chain.size();
    }
  }

  m_surface_intersections = Matrix3Xf(3, sum);
  int i = 0;
  for (auto chain : intersections_chains) {
    auto it = chain.begin();
    SVertex *sv0 = *it;
    m_surface_intersections.col(i++) = sv0->m_pos3D;
    it++;
    while (it != chain.end()) {
      SVertex *sv1 = *it;
      m_surface_intersections.col(i++) = sv1->m_pos3D;
      sv0 = sv1;
      it++;
    }
  }

  for (auto vec : surf_intersect.vector())
    for (auto s : vec)
      delete s;
  remove_edge_property(surf_intersect);
}

void Mesh::projectSVertices(const Matrix4f &model, const Matrix4f &proj,
                            const Vector2i &viewportSize) {
  for (auto v : m_svertices) {
    v->m_pos2D = project(v->m_pos3D, model, proj, viewportSize);
  }
}

/// \returns true if d and e are on the same side of the abc triangle
inline bool Mesh::sameSide(const Eigen::Vector3d &a, const Eigen::Vector3d &b,
                           const Eigen::Vector3d &c, const Eigen::Vector3d &d,
                           const Eigen::Vector3d &e) {
  return (GEO::PCK::orient_3d(a.data(), b.data(), c.data(), d.data()) > 0) ==
         (GEO::PCK::orient_3d(a.data(), b.data(), c.data(), e.data()) > 0);
}

// \returns returns a positive value if d is on the front side of triangle abc
inline double Mesh::frontSide(const Eigen::Vector3d &a,
                              const Eigen::Vector3d &b,
                              const Eigen::Vector3d &c,
                              const Eigen::Vector3d &d) {
  return GEO::PCK::orient_3d(a.data(), b.data(), c.data(), d.data());
}

/// \returns true if c and d are on the same side of the line ab
inline bool Mesh::sameSide(const Eigen::Vector2d &a, const Eigen::Vector2d &b,
                           const Eigen::Vector2d &c, const Eigen::Vector2d &d) {
  return (GEO::PCK::orient_2d(a.data(), b.data(), c.data()) > 0) ==
         (GEO::PCK::orient_2d(a.data(), b.data(), d.data()) > 0);
}

void Mesh::computeContourChains(ContourMode mode) {
  auto is_contour = get_edge_property<real_t>("e:contour");
  auto visited = get_edge_property<real_t>("e:visited");
  if (!visited) {
    visited = add_edge_property<real_t>("e:visited", false);
  } else {
    for (size_t i = 0; i < visited.vector().size(); ++i)
      visited.vector()[i] = false;
  }
  m_chains.clear();
  m_chain_lengths.clear();

  size_t sum = 0;

  for (size_t i = 0; i < m_contour_edges.size(); ++i) {
    Edge e = m_contour_edges[i];
    if (visited[e])
      continue;
    std::list<Edge> chain;
    chain.push_back(e);
    visited[e] = true;
    bool forward = true;
    if (mode == INTERPOLATED_CONTOUR) {
      Halfedge nh = next_halfedge(halfedge(e, 0));
      Edge ne = edge(nh);
      while (ne != e) {
        if (is_contour[ne] >= 0 && !visited[ne]) {
          if (forward)
            chain.push_back(ne);
          else
            chain.push_front(ne);
          nh = next_halfedge(opposite_halfedge(nh));
        } else {
          nh = next_halfedge(nh);
          if (is_boundary(nh)) {
            if (forward) {
              // Chain backward
              forward = false;
              nh = next_halfedge(halfedge(e, 1));
            } else {
              break;
            }
          }
        }
        visited[ne] = true;
        ne = edge(nh);
      }
      if (ne == e)
        chain.push_back(e);
    } else { // MESH CONTOURS
      Halfedge nh = next_halfedge(halfedge(e, 0));
      Edge ne = edge(nh);
      Edge ie = ne;
      while (ne != e) {
        if (is_contour[ne] >= 0 && !visited[ne]) {
          if (forward)
            chain.push_back(ne);
          else
            chain.push_front(ne);
          nh = opposite_halfedge(nh);
          ie = edge(nh);
        } else {
          nh = next_halfedge(opposite_halfedge(nh));
          if (edge(nh) == ie) { // full loop around the one-ring
            if (forward) {
              // Chain backward
              forward = false;
              nh = next_halfedge(halfedge(e, 1));
              ie = edge(nh);
            } else {
              break;
            }
          }
        }
        visited[ne] = true;
        ne = edge(nh);
      }
    }
    m_chains.push_back(chain);
    m_chain_lengths.push_back(chain.size());
    sum += chain.size();
  }

  auto vpositions = get_vertex_property<Vector3f>("v:point");
  auto fedges_e = edge_property<FEdge *>("e:fedge", nullptr);
  auto svertices_e = edge_property<std::vector<SVertex *>>("e:svertex");

  int i = 0;
  if (mode == INTERPOLATED_CONTOUR) {
    auto fedges_f = face_property<std::vector<FEdgeIntersection *>>("f:fedge");
    auto addSvertex = [&](Edge e, const Vector3f &p) {
      SVertex *sv = nullptr;
      if (!svertices_e[e].empty()) {
        for (auto v : svertices_e[e])
          if ((v->m_pos3D - p).squaredNorm() < EPSILON) {
            sv = v;
          }
      }
      if (!sv) {
        sv = new SVertex(p);
        svertices_e[e].push_back(sv);
        m_svertices.push_back(sv);
      }
      if (fedges_e[e]) { // 3D intersection
        sv->m_nature = VertexNature::INTERSECTION_3D | sv->m_nature;
      }
      return sv;
    };
    auto commonFace = [&](Edge e1, Edge e2) {
      Face f1[2] = {face(halfedge(e1, 0)), face(halfedge(e1, 1))};
      Face f2[2] = {face(halfedge(e2, 0)), face(halfedge(e2, 1))};
      for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
          if (f1[i] == f2[j])
            return f1[i];
        }
      }
      std::runtime_error("No face in common");
      return Face();
    };

    m_interpolated_contours = Matrix3Xf(3, sum + m_chains.size() * 2);
    for (auto chain : m_chains) {
      SVertex *sv0, *sv1;
      sv0 = sv1 = nullptr;
      real_t w0;
      Edge e0;
      Vector3f p;
      for (Edge e1 : chain) {
        Vertex i0 = vertex(e1, 0);
        Vertex i1 = vertex(e1, 1);
        real_t w1 = is_contour[e1];
        Vector3f p0 = vpositions[i0];
        Vector3f p1 = vpositions[i1];
        Vector3f p = w1 * p1 + (1.f - w1) * p0;
        m_interpolated_contours.col(i++) = p;
        if (!sv0) { // first vertex
          sv0 = addSvertex(e1, p);
        } else {
          sv1 = addSvertex(e1, p);

          // test if 3D intersection
          bool intersection_found = false;
          Face f = commonFace(e0, e1);
          for (FEdgeIntersection *fe : fedges_f[f]) {
            Vector3f s0 = p0, s1 = p1, s2;
            Eigen::Vector2d pC, pD;
            if (vertex(e0, 0) == i0) {
              s2 = vpositions[vertex(e0, 1)];
              pC = Eigen::Vector2d(w1, 0);
              pD = Eigen::Vector2d(0, w0);
            } else if (vertex(e0, 0) == i1) {
              s2 = vpositions[vertex(e0, 1)];
              std::swap(s0, s1);
              pC = Eigen::Vector2d(1.0 - w1, 0);
              pD = Eigen::Vector2d(0, w0);
            } else if (vertex(e0, 1) == i0) {
              s2 = vpositions[vertex(e0, 0)];
              pC = Eigen::Vector2d(w1, 0);
              pD = Eigen::Vector2d(0, 1.0 - w0);
            } else if (vertex(e0, 1) == i1) {
              s2 = vpositions[vertex(e0, 0)];
              std::swap(s0, s1);
              pC = Eigen::Vector2d(1.0 - w1, 0);
              pD = Eigen::Vector2d(0, 1.0 - w0);
            } else {
              std::runtime_error("Invalid contour chain");
            }
            Vector3f v0 = s1 - s0, v1 = s2 - s0,
                     v2 = fe->m_vertexA->m_pos3D - s0;
            double d00 = v0.squaredNorm();
            double d01 = v0.dot(v1);
            double d11 = v1.squaredNorm();
            double d20 = v2.dot(v0);
            double d21 = v2.dot(v1);
            double inv_denom = 1.0 / (d00 * d11 - d01 * d01);
            Eigen::Vector2d pA =
                inv_denom *
                Eigen::Vector2d(d11 * d20 - d01 * d21, d00 * d21 - d01 * d20);
            v2 = fe->m_vertexB->m_pos3D - s0;
            d20 = v2.dot(v0);
            d21 = v2.dot(v1);
            Eigen::Vector2d pB =
                inv_denom *
                Eigen::Vector2d(d11 * d20 - d01 * d21, d00 * d21 - d01 * d20);
            double s, t;
            if (intersect2dSeg2dSegParametric(pA.data(), pB.data(), pC.data(),
                                              pD.data(), s, t, EPSILON) == 1) {
              Eigen::Vector2d uv = pB * s + pA * (1.0 - s);

              // Split fegdes at intersection
              SVertex *sv_intersec =
                  new SVertex(v0 * uv.x() + v1 * uv.y() + s0);
              sv_intersec->m_nature = VertexNature::INTERSECTION_3D;
              m_svertices.push_back(sv_intersec);
              FEdgeIntersection *f1, *f2, *f3, *f4;
              fe->m_vertexA->m_fedges.remove(fe);
              fe->m_vertexB->m_fedges.remove(fe);
              f1 = new FEdgeIntersection(fe->m_vertexA, sv_intersec, fe->m_f1,
                                         fe->m_f2);
              fe->m_vertexA->m_fedges.push_back(f1);
              sv_intersec->m_fedges.push_back(f1);
              f2 = new FEdgeIntersection(sv_intersec, fe->m_vertexB, fe->m_f1,
                                         fe->m_f2);
              fe->m_vertexB->m_fedges.push_back(f2);
              sv_intersec->m_fedges.push_back(f2);
              f3 = new FEdgeIntersection(sv0, sv_intersec, fe->m_f1, fe->m_f2);
              sv0->m_fedges.push_back(f3);
              sv_intersec->m_fedges.push_back(f3);
              f4 = new FEdgeIntersection(sv_intersec, sv1, fe->m_f1, fe->m_f2);
              sv1->m_fedges.push_back(f4);
              sv_intersec->m_fedges.push_back(f4);
              auto it = std::find(m_fedges.begin(), m_fedges.end(), fe);
              if (it != m_fedges.end())
                m_fedges.erase(it);
              m_fedges.push_back(f1);
              m_fedges.push_back(f2);
              m_fedges.push_back(f3);
              m_fedges.push_back(f4);
              intersection_found = true;
            }
          }
          if (!intersection_found) {
            FEdge *fedge = new FEdgeSmooth(sv0, sv1, e0, w0, e1, w1,
                                           EdgeNature::SMOOTH_CONTOUR);
            m_fedges.push_back(fedge);
            sv0->m_fedges.push_back(fedge);
            sv1->m_fedges.push_back(fedge);
          }
          sv0 = sv1;
        }
        w0 = w1;
        e0 = e1;
      }
    }
  } else { // MESH CONTOURS
    auto svertices_v = vertex_property<SVertex *>("v:svertex", nullptr);
    auto concave = get_edge_property<bool>("e:concave");
    // create a fedge between two vertices of the mesh
    auto add_fedge = [&](Edge edge, int i0, SVertex *sv0) {
      if (!sv0) {
        Vertex v0 = vertex(edge, i0);
        if (svertices_v[v0]) {
          sv0 = svertices_v[v0];
          // 3D intersection
          sv0->m_nature = VertexNature::INTERSECTION_3D | sv0->m_nature;
        } else {
          sv0 = new SVertex(vpositions[v0]);
          svertices_v[v0] = sv0;
          m_svertices.push_back(sv0);
        }
      }
      Vertex v1 = vertex(edge, (i0 + 1) % 2);
      SVertex *sv1;
      if (svertices_v[v1]) {
        sv1 = svertices_v[v1];
        if (sv1->m_fedges.size() >= 2) { // 3D intersection
          sv1->m_nature = VertexNature::INTERSECTION_3D | sv1->m_nature;
        }
      } else {
        sv1 = new SVertex(vpositions[v1]);
        svertices_v[v1] = sv1;
        m_svertices.push_back(sv1);
      }
      FEdge *fedge = new FEdgeSharp(sv0, sv1, edge, EdgeNature::SHARP_CONTOUR);
      m_fedges.push_back(fedge);
      fedges_e[edge] = fedge;
      sv0->m_fedges.push_back(fedge);
      sv1->m_fedges.push_back(fedge);

      if (!svertices_e[edge].empty()) { // 3D intersection
        for (auto v : svertices_e[edge]) {
          v->m_nature = VertexNature::INTERSECTION_3D | v->m_nature;
        }
      }
      return sv1;
    };

    m_chain_indices.resize(sum + m_chains.size() * 3);
    for (auto chain : m_chains) {
      Edge e = chain.front();
      SVertex *sv = nullptr;
      if (chain.size() == 1) {
        m_chain_indices[i++] = vertex(e, 0).idx();
        m_chain_indices[i++] = vertex(e, 1).idx();
        add_fedge(e, 0, sv);
      } else {
        auto it = chain.begin();
        it++;
        Edge ne = (*it);
        Vertex prev;
        bool prevConcave;
        if (vertex(e, 0) == vertex(ne, 0) || vertex(e, 0) == vertex(ne, 1)) {
          m_chain_indices[i++] = vertex(e, 1).idx();
          prev = vertex(e, 0);
          sv = add_fedge(e, 1, sv);
        } else if (vertex(e, 1) == vertex(ne, 0) ||
                   vertex(e, 1) == vertex(ne, 1)) {
          m_chain_indices[i++] = vertex(e, 0).idx();
          prev = vertex(e, 1);
          sv = add_fedge(e, 0, sv);
        } else {
          std::runtime_error("Invalid contour chain");
        }
        prevConcave = concave[e];
        m_chain_indices[i++] = prev.idx();
        while (it != chain.end()) {
          e = *it;
          if (prevConcave != concave[e]) { // curtain fold
            sv->m_nature = VertexNature::CONTOUR_CURTAIN_FOLD | sv->m_nature;
            m_debug_points.push_back(sv->m_pos3D);
            m_point_colors.push_back(curtainFoldColor.head<3>().cast<real_t>());
          }
          prevConcave = concave[e];
          if (prev != vertex(e, 0)) {
            prev = vertex(e, 0);
            sv = add_fedge(e, 1, sv);
          } else {
            prev = vertex(e, 1);
            sv = add_fedge(e, 0, sv);
          }
          m_chain_indices[i++] = prev.idx();
          it++;
        }
        // m_chain_indices[i++] = prev.idx();
      }
    }
  }
}

void Mesh::computeSweepLineIntersections(ContourMode mode,
                                         const Matrix4f &projectionMatrix,
                                         const Matrix4f &viewMatrix,
                                         const Vector2i &viewport) {
  typedef Segment<FEdge *> segment;
  typedef Intersection<segment> intersection;

  // we only want image-space intersections between any pair of edges except
  // for smooth-sharp pairs that intersect on the surface
  struct silhouette_binary_rule_no_same_face
      : public binary_rule<segment, segment> {
    silhouette_binary_rule_no_same_face() : binary_rule<segment, segment>() {}
    virtual bool operator()(segment &s1, segment &s2) {
      FEdge *e1 = s1.edge();
      FEdge *e2 = s2.edge();

      if ((e1->isSmooth() && e2->isSmooth()) ||
          (!e1->isSmooth() && !e2->isSmooth()))
        return true;

      if (e1->isSmooth()) {
        if (e1->m_vertexA->m_nature & VertexNature::INTERSECTION_3D ||
            e1->m_vertexB->m_nature & VertexNature::INTERSECTION_3D)
          return false;
      } else {
        if (e2->m_vertexA->m_nature & VertexNature::INTERSECTION_3D ||
            e2->m_vertexB->m_nature & VertexNature::INTERSECTION_3D)
          return false;
      }
      return true;
    }
  };

  vector<segment *> segments;

  sort(m_svertices.begin(), m_svertices.end(), VertexCompare(EPSILON));

  for (auto fedge : m_fedges) {
    segments.push_back(new segment(fedge, fedge->m_vertexA, fedge->m_vertexB));
    fedge->m_segment = segments.back();
  }

  SweepLine<FEdge *> SL;
  vector<segment *> vsegments;
  silhouette_binary_rule_no_same_face sbr;
  // Iterate over every sorted curve vertices
  for (auto v : m_svertices) {
    // Add new fedges
    for (auto fedge : v->m_fedges) {
      assert(fedge);
      vsegments.push_back(fedge->m_segment);
    }
    SL.process(v->m_pos2D.cast<double>().head<2>(), vsegments, sbr);
    vsegments.clear();
  }

  // retrieve the intersections
  vector<intersection *> &intersections = SL.intersections();
  m_2D_intersections.second = intersections.size();

  auto addIntersectionSVertex = [&](FEdge *e, real_t t) {
    SVertex *svA = e->m_vertexA;
    SVertex *svB = e->m_vertexB;
    Vector3f pos2D = svA->m_pos2D + t * (svB->m_pos2D - svA->m_pos2D);
    Vector3f pos3D = unproject(pos2D, viewMatrix, projectionMatrix, viewport);
    // create new SVertex
    SVertex *sv_intersec = new SVertex(pos3D);
    sv_intersec->m_pos2D = pos2D;
    sv_intersec->m_nature =
        sv_intersec->m_nature | VertexNature::INTERSECTION_2D;
    m_svertices.push_back(sv_intersec);
    return sv_intersec;
  };

  // create new svertices
  for (auto intersect : intersections) {
    intersect->m_vertexA =
        addIntersectionSVertex(intersect->m_edgeA->edge(), intersect->m_tA);
    intersect->m_vertexB =
        addIntersectionSVertex(intersect->m_edgeB->edge(), intersect->m_tB);

    // create an FEdge between them
    FEdge *new_edge = new FEdge(intersect->m_vertexA, intersect->m_vertexB,
                                EdgeNature::IMAGE_INTERSECTION);
    intersect->m_vertexA->m_fedges.push_back(new_edge);
    intersect->m_vertexB->m_fedges.push_back(new_edge);
    m_fedges.push_back(new_edge);
  }

  auto addFEdge = [&](FEdge *e, SVertex *prev_sv, SVertex *curr_sv) {
    FEdge *new_edge;
    if (e->m_nature & EdgeNature::SHARP_CONTOUR ||
        e->m_nature & EdgeNature::BOUNDARY) {
      FEdgeSharp *fe_s = dynamic_cast<FEdgeSharp *>(e);
      new_edge = new FEdgeSharp(prev_sv, curr_sv, fe_s->m_edge, e->m_nature);
    } else if (e->m_nature & EdgeNature::SMOOTH_CONTOUR) {
      FEdgeSmooth *fe_s = dynamic_cast<FEdgeSmooth *>(e);
      new_edge = new FEdgeSmooth(prev_sv, curr_sv, fe_s->m_edgeA, fe_s->m_tA,
                                 fe_s->m_edgeB, fe_s->m_tB, e->m_nature);
    } else if (e->m_nature & EdgeNature::SURFACE_INTERSECTION) {
      FEdgeIntersection *fe_i = dynamic_cast<FEdgeIntersection *>(e);
      new_edge =
          new FEdgeIntersection(prev_sv, curr_sv, fe_i->m_f1, fe_i->m_f2);
    }
    curr_sv->m_fedges.push_back(new_edge);
    m_fedges.push_back(new_edge);
  };

  // retrieve the intersected edges
  std::vector<segment *> &iedges = SL.intersectedEdges();
  // split egdes at intersection
  for (auto s : iedges) {
    vector<intersection *> &eIntersections = s->intersections();
    // first need to sort these intersections along the edge
    std::sort(eIntersections.begin(), eIntersections.end(),
              less_Intersection<FEdge>(s));

    FEdge *e = s->edge();
    // retrieve the new svertices for all intersections on this edge
    vector<SVertex *> edgeSVertices;
    for (auto intersect : eIntersections) {
      if (e == intersect->m_edgeA->edge()) {
        edgeSVertices.push_back(intersect->m_vertexA);
      } else {
        assert(e == intersect->m_edgeB->edge());
        edgeSVertices.push_back(intersect->m_vertexB);
      }
    }

    // remove old FEdge from its extremities
    SVertex *svA = e->m_vertexA;
    SVertex *svB = e->m_vertexB;
    svA->m_fedges.remove(e);
    svB->m_fedges.remove(e);
    // chain new vertices
    SVertex *prev_sv = svA;
    for (auto curr_sv : edgeSVertices) {
      addFEdge(e, prev_sv, curr_sv);
      prev_sv = curr_sv;
    }
    addFEdge(e, prev_sv, svB);

    // delete the old FEdge 
    auto it = std::find(m_fedges.begin(), m_fedges.end(), e);
    if (it != m_fedges.end())
      m_fedges.erase(it);
    delete e;
  }

  // clean-up
  for (auto fedge : m_fedges)
    fedge->m_segment = nullptr;
  for (auto s : segments)
    delete s;
}

void Mesh::computeContourVisibility(const Vector3f &view_pos,
                                    ContourMode mode) {

  if (!m_bvh) {
    m_bvh = new BVH(
        &m_indices, &get_vertex_property<Vector3f>("v:point").vector(),
        &get_vertex_property<Vector3f>("v:normal").vector(), boundingBox());
    m_bvh->build();
  }

  auto concave = get_edge_property<bool>("e:concave");
  ray_intersections.clear();

  if (mode != INTERPOLATED_CONTOUR) {
    std::vector<Vector2u> invisible_indices;
    invisible_indices.reserve(0.5 * m_mesh_contours_indices.size());
    std::vector<Vector2u> visible_indices;
    visible_indices.reserve(0.5 * m_mesh_contours_indices.size());

    auto visibility = get_edge_property<bool>("e:visible");
    if (!visibility)
      visibility = add_edge_property<bool>("e:visible");

    auto vpositions = get_vertex_property<Vector3f>("v:point");
    for (uint32_t i = 0; i < m_contour_edges.size(); ++i) {
      Edge e = m_contour_edges[i];
      const Vector3f &a = vpositions[vertex(e, 0)];
      const Vector3f &b = vpositions[vertex(e, 1)];

      if (concave[e]) {
        //  ray_intersections.push_back(0.5f * (a + b));
        continue;
      }

      Vector3f dir = real_t(0.5) * (a + b) - view_pos;
      uint32_t idx;
      real_t t;
      visibility[e] = true;
      Ray ray = Ray(view_pos, dir.normalized(), 0, dir.norm());
      bool hit = m_bvh->rayIntersect(ray, idx, t);
      if (hit && int(idx) != face(e, 0).idx() && int(idx) != face(e, 1).idx()) {
        ray_intersections.push_back(ray(t));
        visibility[e] = false;
        invisible_indices.push_back(
            Vector2u(vertex(e, 0).idx(), vertex(e, 1).idx()));
      } else {
        visible_indices.push_back(
            Vector2u(vertex(e, 0).idx(), vertex(e, 1).idx()));
      }
    }

    m_mesh_contours_indices =
        MatrixXu(2, visible_indices.size() + invisible_indices.size());
    for (size_t i = 0; i < visible_indices.size(); ++i)
      m_mesh_contours_indices.col(i) = visible_indices.at(i);
    for (size_t i = 0; i < invisible_indices.size(); ++i)
      m_mesh_contours_indices.col(i + visible_indices.size()) =
          invisible_indices.at(i);
    m_visible = 2 * visible_indices.size();
  } else { // INTERPOLATED CONTOURS
    std::vector<std::pair<Vector3f, Vector3f>> invisible_edges;
    invisible_edges.reserve(0.25 * m_interpolated_contours.cols());
    std::vector<std::pair<Vector3f, Vector3f>> visible_edges;
    visible_edges.reserve(0.25 * m_interpolated_contours.cols());
    for (uint32_t i = 0; i < m_interpolated_contours.cols() / 2; ++i) {
      const Vector3f &p1 = m_interpolated_contours.col(2 * i);
      const Vector3f &p2 = m_interpolated_contours.col(2 * i + 1);
      Vector3f dir = 0.5f * (p1 + p2) - view_pos;
      uint32_t idx;
      real_t t;
      Ray ray = Ray(view_pos, dir.normalized(), 0, dir.norm());
      bool hit = m_bvh->rayIntersect(ray, idx, t);
      if (hit && int(idx) != m_interpolated_contour_faces[i]) {
        ray_intersections.push_back(ray(t));
        invisible_edges.push_back(std::make_pair(p1, p2));
      } else {
        visible_edges.push_back(std::make_pair(p1, p2));
      }
    }
    for (size_t i = 0; i < visible_edges.size(); ++i) {
      m_interpolated_contours.col(i * 2) = visible_edges.at(i).first;
      m_interpolated_contours.col(i * 2 + 1) = visible_edges.at(i).second;
    }
    for (size_t i = 0; i < invisible_edges.size(); ++i) {
      m_interpolated_contours.col((i + visible_edges.size()) * 2) =
          invisible_edges.at(i).first;
      m_interpolated_contours.col((i + visible_edges.size()) * 2 + 1) =
          invisible_edges.at(i).second;
    }
    m_visible = 2 * visible_edges.size();
  }
}

Matrix3Xf Mesh::get_contours(ContourMode mode, bool visible_only) {
  if (mode == INTERPOLATED_CONTOUR) {
    if (visible_only)
      return m_interpolated_contours.topLeftCorner(3, m_visible);
    return get_interpolated_contours();
  } else {
    auto vpositions = get_vertex_property<Vector3f>("v:point");

    Matrix3Xf contours = Matrix3Xf(3, m_chain_indices.size());
    for (size_t i = 0; i < m_chain_indices.size(); ++i) {
      Vertex v = Vertex(m_chain_indices[i]);
      contours.col(i) = vpositions[v];
    }
    return contours;
  }
  return Matrix3Xf();
}

Matrix3Xf Mesh::get_boundaries(bool visible_only) {
  auto vpositions = get_vertex_property<Vector3f>("v:point");

  Matrix3Xf boundaries = Matrix3Xf(3, m_boundary_indices.size());
  for (size_t i = 0; i < m_boundary_indices.size(); ++i) {
    Vertex v = Vertex(m_boundary_indices[i]);
    boundaries.col(i) = vpositions[v];
  }

  return boundaries;
}

std::vector<Vector3f> &Mesh::get_debug_points(nature_t nature) {
  using namespace VertexNature;

  m_debug_points.clear();
  m_point_colors.clear();
  m_2D_intersections.first = -1;

  for (auto sv : m_svertices) {
    if (!sv)
      continue;
    nature_t sv_nature = sv->m_nature;
    if (nature & INTERSECTION_3D && sv_nature & INTERSECTION_3D) {
      m_debug_points.push_back(sv->m_pos3D);
      m_point_colors.push_back(intersection3DColor.head<3>().cast<real_t>());
    }
    if (nature & INTERSECTION_2D && sv_nature & INTERSECTION_2D) {
      if (m_2D_intersections.first == -1)
        m_2D_intersections.first = m_debug_points.size();
      m_debug_points.push_back(sv->m_pos3D);
      m_point_colors.push_back(intersection2DColor.head<3>().cast<real_t>());
    }
    if (nature & BOUNDARY_CURTAIN_FOLD && sv_nature & BOUNDARY_CURTAIN_FOLD) {
      m_debug_points.push_back(sv->m_pos3D);
      m_point_colors.push_back(curtainFoldColor.head<3>().cast<real_t>());
    }
    if (nature & CONTOUR_CURTAIN_FOLD && sv_nature & CONTOUR_CURTAIN_FOLD) {
      m_debug_points.push_back(sv->m_pos3D);
      m_point_colors.push_back(curtainFoldColor.head<3>().cast<real_t>());
    }
    if (nature & BIFURCATION && sv_nature & BIFURCATION) {
      m_debug_points.push_back(sv->m_pos3D);
      m_point_colors.push_back(bifurcationColor.head<3>().cast<real_t>());
    }
  }

  m_debug_points.reserve(m_debug_points.size() + ray_intersections.size());
  m_point_colors.reserve(m_point_colors.size() + ray_intersections.size());
  for (size_t i = 0; i < ray_intersections.size(); i++) {
    m_debug_points.push_back(ray_intersections.at(i));
    m_point_colors.push_back(Vector3f(0.9, 0.0, 0.0));
  }

  return m_debug_points;
}
