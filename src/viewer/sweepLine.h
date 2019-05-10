#ifndef SWEEPLINE_H
#define SWEEPLINE_H

#include "common.h"
#include "mesh.h"
#include "viewgraph.h"

using namespace std;

/*! Class to define the intersection between two segments*/
template <class Edge> class Intersection {
public:
  template <class EdgeClass>
  Intersection(EdgeClass *eA, double ta, EdgeClass *eB, double tb) {
    m_edgeA = eA;
    m_edgeB = eB;
    m_tA = ta;
    m_tB = tb;
  }

  /*! returns the parameter giving the intersection, for the edge iEdge */
  double getParameter(Edge *iEdge) {
    if (iEdge == m_edgeA)
      return m_tA;
    if (iEdge == m_edgeB)
      return m_tB;
    return 0;
  }

public:
  Edge *m_edgeA; // first segment
  Edge *m_edgeB; // second segment
  double m_tA; // parameter defining the intersection point with respect to the
               // segment EdgeA.
  double m_tB; // parameter defining the intersection point with respect to the
               // segment EdgeB.

  SVertex *m_vertexA, *m_vertexB; // the new vertices
};

template <class Edge> class Segment {
  bool less(const Eigen::Vector2d &v1, const Eigen::Vector2d &v2) {
    for (unsigned int i = 0; i < 2; i++) {
      if (v1[i] < v2[i])
        return true;
      if (v1[i] > v2[i])
        return false;
    }
    return false;
  }

public:
  Segment(Edge &s, SVertex *sA, SVertex *sB) {
    m_edge = s;
    m_sA = sA;
    m_sB = sB;
    Eigen::Vector2d iA = m_sA->m_pos2D.cast<double>().head<2>();
    Eigen::Vector2d iB = m_sB->m_pos2D.cast<double>().head<2>();
    if (less(iA, iB)) {
      m_pA = iA;
      m_pB = iB;
      m_sameOrder = true;
    } else {
      m_pA = iB;
      m_pB = iA;
      m_sameOrder = false;
    }
  }

  ~Segment() { m_intersections.clear(); }

  inline Eigen::Vector2d operator[](const unsigned short int &i) const {
    return i % 2 == 0 ? m_pA : m_pB;
  }

  inline bool operator==(const Segment<Edge> &iBrother) {
    return (m_edge == iBrother.m_edge);
  }

  /* Adds an intersection for this segment */
  inline void addIntersection(Intersection<Segment<Edge>> *i) {
    m_intersections.push_back(i);
  }

  /*! Checks for a common vertex with another edge */
  inline bool commonVertex(const Segment<Edge> &s) {
    if (m_sA == s.m_sA || m_sA == s.m_sB)
      return true;
    if (m_sB == s.m_sA || m_sB == s.m_sB)
      return true;
    return false;
  }

  inline vector<Intersection<Segment<Edge>> *> &intersections() {
    return m_intersections;
  }
  inline bool sameOrder() { return m_sameOrder; }
  inline Edge &edge() { return m_edge; }

private:
  Edge m_edge;
  SVertex *m_sA, *m_sB;
  Eigen::Vector2d m_pA;
  Eigen::Vector2d m_pB;
  // list of intersections parameters
  std::vector<Intersection<Segment<Edge>> *> m_intersections;
  bool m_sameOrder; // true if A and B are in the same order than m_edge.A and
                    // m_edge.B. false otherwise.
};

template <class Edge>
struct less_Intersection
    : public binary_function<Intersection<Segment<Edge *>> *,
                             Intersection<Segment<Edge *>> *, bool> {
  Segment<Edge *> *edge;

  less_Intersection(Segment<Edge *> *iEdge)
      : binary_function<Intersection<Segment<Edge *>> *,
                        Intersection<Segment<Edge *>> *, bool>() {
    edge = iEdge;
  }

  bool operator()(Intersection<Segment<Edge *>> *x,
                  Intersection<Segment<Edge *>> *y) {
    double tx = x->getParameter(edge);
    double ty = y->getParameter(edge);
    if (tx > ty)
      return true;
    return false;
  }
};

/*! defines a binary function that can be overload
 *  by the user to specify at each condition
 *  the intersection between 2 edges must be computed
 */
template <class T1, class T2> struct binary_rule {
  binary_rule() {}
  template <class T3, class T4>
  binary_rule(const binary_rule<T3, T4> &brother) {}
  virtual ~binary_rule() {}

  virtual bool operator()(T1 &, T2 &) { return true; }
};

extern int intersect2dSeg2dSegParametric(double p1[2], double p2[2],
                                         double p3[2], double p4[2], double &t,
                                         double &u, double epsilon);

template <class Edge> class SweepLine {

public:
  SweepLine() {}

  ~SweepLine() {
    for (auto i = m_intersections.begin(), iend = m_intersections.end();
         i != iend; i++) {
      delete (*i);
    }
    m_intersections.clear();
    m_intersectedEdges.clear();

    m_set.clear();
  }

  // it is assumed that the points have been sorted according to X coordinate
  // This function is called for each point, in order. "segments" are all
  // segments that touch this point. Hence, if a segment ends at p, then it can
  // be removed from the active set of segments.

  inline void process(const Eigen::Vector2d &p,
                      const vector<Segment<Edge> *> &segments,
                      binary_rule<Segment<Edge>, Segment<Edge>> &binrule) {
    // first we remove the segments that need to be removed and then
    // we add the segments to add
    vector<Segment<Edge> *> toadd;
    typename vector<Segment<Edge> *>::iterator s, send;
    for (auto s : segments) {
      if (p == (*s)[0]) {
        toadd.push_back(s);
      } else {
        remove(s);
      }
    }
    for (auto s : toadd) {
      add(s, binrule);
    }
  }

  inline void add(Segment<Edge> *S,
                  binary_rule<Segment<Edge>, Segment<Edge>> &binrule) {
    double t, u;
    Eigen::Vector2d v0, v1, v2, v3;
    if (S->sameOrder()) {
      v0 = ((*S)[0]);
      v1 = ((*S)[1]);
    } else {
      v0 = ((*S)[1]);
      v1 = ((*S)[0]);
    }
    for (auto currentS : m_set) {
      if (!binrule(*S, *currentS))
        continue;

      if (S->commonVertex(*currentS))
        continue; // the two edges have a common vertex => no need to check

      if (currentS->sameOrder()) {
        v2 = ((*currentS)[0]);
        v3 = ((*currentS)[1]);
      } else {
        v2 = ((*currentS)[1]);
        v3 = ((*currentS)[0]);
      }

      if (intersect2dSeg2dSegParametric(v0.data(), v1.data(), v2.data(),
                                        v3.data(), t, u, EPSILON) == 1) {
        // create the intersection
        Intersection<Segment<Edge>> *inter =
            new Intersection<Segment<Edge>>(S, t, currentS, u);
        // assert(t > EPSILON && u > EPSILON);
        // add it to the intersections list
        m_intersections.push_back(inter);
        // add this intersection to the first edge intersections list
        S->addIntersection(inter);
        // add this intersection to the second edge intersections list
        currentS->addIntersection(inter);
      }
    }
    // add the added segment to the list of active segments
    m_set.push_back(S);
  }

  inline void remove(Segment<Edge> *s) {
    if (s->intersections().size() > 0)
      m_intersectedEdges.push_back(s);
    m_set.remove(s);
  }

  vector<Segment<Edge> *> &intersectedEdges() { return m_intersectedEdges; }

  vector<Intersection<Segment<Edge>> *> &intersections() {
    return m_intersections;
  }

private:
  // set of active edges for a given position of the sweep line
  std::list<Segment<Edge> *> m_set;
  // the list of intersected edges
  std::vector<Segment<Edge> *> m_intersectedEdges;
  // the list of all intersections.
  std::vector<Intersection<Segment<Edge>> *> m_intersections;
};

#endif