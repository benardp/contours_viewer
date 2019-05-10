/*
    bvh.h -- bounding volume hierarchy for fast ray-intersection queries

    This file is part of the implementation of

        Instant Field-Aligned Meshes
        Wenzel Jakob, Daniele Panozzo, Marco Tarini, and Olga Sorkine-Hornung
        In ACM Transactions on Graphics (Proc. SIGGRAPH Asia 2015)

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE.txt file.
*/

#ifndef BVH_H
#define BVH_H

#include "aabb.h"

// #if defined(SINGLE_PRECISION)
  typedef uint32_t uintp_t;
// #else
//   typedef uint64_t uintp_t;
// #endif

/* BVH node in 32 bytes */
struct BVHNode {
  union {
    struct {
      unsigned flag : 1;
      uintp_t size : 31;
      uintp_t start;
    } leaf;

    struct {
      uintp_t unused;
      uintp_t rightChild;
    } inner;
  };
  AABB aabb;

  inline bool isLeaf() const { return leaf.flag == 1; }

  inline bool isInner() const { return leaf.flag == 0; }

  inline bool isUnused() const {
    return inner.unused == 0 && inner.rightChild == 0;
  }

  inline uintp_t start() const { return leaf.start; }

  inline uintp_t end() const { return leaf.start + leaf.size; }
};

class BVH {
  friend struct BVHBuildTask;
  /* Cost values for BVH surface area heuristic */
  enum { T_aabb = 1, T_tri = 1 };

public:
  BVH(std::vector<uint32_t> *F, std::vector<Vector3f> *V,
      std::vector<Vector3f> *N, const AABB &aabb);

  ~BVH();

  void build();

  void printStatistics() const;

  bool rayIntersect(Ray ray) const;

  bool rayIntersect(Ray ray, uint32_t &idx, real_t &t,
                    Vector2f *uv = nullptr) const;

  bool rayIntersect(Ray ray, std::vector<uint32_t> &idx,
                    std::vector<Vector2f> *uvs = nullptr) const;

  void findNearestWithRadius(const Vector3f &p, real_t radius,
                             std::vector<uint32_t> &result,
                             bool includeSelf = false) const;

  uint32_t findNearest(const Vector3f &p, real_t &radius,
                       bool includeSelf = false) const;

  void findKNearest(const Vector3f &p, uint32_t k, real_t &radius,
                    std::vector<std::pair<real_t, uint32_t>> &result,
                    bool includeSelf = false) const;

  void findKNearest(const Vector3f &p, const Vector3f &N, uint32_t k,
                    real_t &radius,
                    std::vector<std::pair<real_t, uint32_t>> &result,
                    real_t angleThresh = 30, bool includeSelf = false) const;

protected:
  bool rayIntersectTri(const Ray &ray, uint32_t i, real_t &t,
                       Vector2f &uv) const;
  std::pair<real_t, uintp_t> statistics(uintp_t node_idx = 0) const;

protected:
  std::vector<BVHNode> mNodes;
  std::vector<uintp_t> mIndices;
  std::vector<uint32_t> *mF;
  std::vector<Vector3f> *mV;
  std::vector<Vector3f> *mN;
};

#endif