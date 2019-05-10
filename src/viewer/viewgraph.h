#ifndef VIEWGRAPH_H
#define VIEWGRAPH_H

#include "common.h"

#include <surface_mesh.h>

using Edge = surface_mesh::Surface_mesh::Edge;
using Face = surface_mesh::Surface_mesh::Face;
using Vertex = surface_mesh::Surface_mesh::Vertex;

struct FEdge;

struct SVertex {
  SVertex(const Vector3f &pos3D, nature_t nature = VertexNature::S_VERTEX)
      : m_pos3D(pos3D), m_nature(nature) {}

  Vector3f m_pos3D, m_pos2D;
  nature_t m_nature;
  std::list<FEdge *> m_fedges;
};

template <class Edge> class Segment;

struct FEdge {
  FEdge(SVertex *vA, SVertex *vB, nature_t n)
      : m_vertexA(vA), m_vertexB(vB), m_nature(n), m_segment(nullptr) {}
  SVertex *m_vertexA, *m_vertexB;

  virtual ~FEdge() {}
  virtual bool isSmooth() const { return false; };

  nature_t m_nature;
  Segment<FEdge *> *m_segment;
};

struct FEdgeSharp : public FEdge {
  FEdgeSharp(SVertex *vA, SVertex *vB, Edge e, nature_t n)
      : FEdge(vA, vB, n), m_edge(e) {}

  virtual ~FEdgeSharp() {}
  virtual bool isSmooth() const { return false; }

  Edge m_edge;
};

struct FEdgeSmooth : public FEdge {
  FEdgeSmooth(SVertex *vA, SVertex *vB, Edge eA, real_t tA, Edge eB, real_t tB,
              nature_t n)
      : FEdge(vA, vB, n), m_edgeA(eA), m_edgeB(eB), m_tA(tA), m_tB(tB) {}

  virtual ~FEdgeSmooth() {}
  virtual bool isSmooth() const { return true; }

  Edge m_edgeA, m_edgeB;
  real_t m_tA, m_tB;
};

struct FEdgeIntersection : public FEdge {
  FEdgeIntersection(SVertex *vA, SVertex *vB, Face f1, Face f2)
      : FEdge(vA, vB, EdgeNature::SURFACE_INTERSECTION), m_f1(f1), m_f2(f2) {}

  virtual ~FEdgeIntersection() {}
  virtual bool isSmooth() const { return true; }

  Face m_f1, m_f2;
};

class VertexCompare {
public:
  VertexCompare(real_t eps) : m_epsilon(eps) {}

  bool operator()(SVertex *x, SVertex *y) const {
    const Vector3f &A = x->m_pos2D;
    const Vector3f &B = y->m_pos2D;
    for (unsigned int i = 0; i < 3; i++) {
      if ((fabs(A[i] - B[i])) < m_epsilon)
        continue;
      if (A[i] < B[i])
        return true;
      if (A[i] > B[i])
        return false;
    }
    return false;
  }

private:
  real_t m_epsilon;
};

#endif
