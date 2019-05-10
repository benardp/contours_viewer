/*
    bvh.cpp -- bounding volume hierarchy for fast ray-intersection queries

    This file is part of the implementation of

        Instant Field-Aligned Meshes
        Wenzel Jakob, Daniele Panozzo, Marco Tarini, and Olga Sorkine-Hornung
        In ACM Transactions on Graphics (Proc. SIGGRAPH Asia 2015)

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE.txt file.
*/

#include "bvh.h"

using namespace std;

struct Bins {
  static const int BIN_COUNT = 8;
  Bins() { memset(counts, 0, sizeof(uintp_t) * BIN_COUNT); }
  uintp_t counts[BIN_COUNT];
  AABB bounds[BIN_COUNT];
};

struct BVHBuildTask {
  enum { SERIAL_THRESHOLD = 32 };
  BVH &bvh;
  uintp_t node_idx;
  uintp_t *start, *end;
  uintp_t *temp;

  /**
   * Create a new build task
   *
   * \param bvh
   *    Reference to the underlying BVH
   *
   * \param node_idx
   *    Index of the BVH node that should be built
   *
   * \param start
   *    Start pointer into a list of triangle indices to be processed
   *
   * \param end
   *    End pointer into a list of triangle indices to be processed
   *
   *  \param temp
   *    Pointer into a temporary memory region that can be used for
   *    construction purposes. The usable length is <tt>end-start</tt>
   *    unsigned integers.
   */
  BVHBuildTask(BVH &bvh, uintp_t node_idx, uintp_t *start, uintp_t *end,
               uintp_t *temp)
      : bvh(bvh), node_idx(node_idx), start(start), end(end), temp(temp) {}

  static void execute_serially(BVH &bvh, uintp_t node_idx, uintp_t *start,
                               uintp_t *end, real_t *temp) {
    uintp_t size = end - start;
    BVHNode &node = bvh.mNodes[node_idx];
    const MapXu F = MapXu(bvh.mF->data(), 3, bvh.mF->size());
    const std::vector<Vector3f> &V = *bvh.mV;
    real_t best_cost = BVH::T_tri * size;
    int64_t best_index = -1, best_axis = -1;
    real_t *left_areas = (real_t *)temp;

    for (int axis = 0; axis < 3; ++axis) {

      std::sort(start, end, [&](uintp_t f1, uintp_t f2) {
        return ((V[F(0, f1)](axis) + V[F(1, f1)](axis) + V[F(2, f1)](axis)) <
                (V[F(0, f2)](axis) + V[F(1, f2)](axis) + V[F(2, f2)](axis)));
      });

      AABB aabb;
      for (uintp_t i = 0; i < size; ++i) {
        uintp_t f = *(start + i);
        aabb.expandBy(V[F(0, f)]);
        aabb.expandBy(V[F(1, f)]);
        aabb.expandBy(V[F(2, f)]);
        left_areas[i] = (real_t)aabb.surfaceArea();
      }
      if (axis == 0)
        node.aabb = aabb;

      aabb.clear();

      real_t tri_factor = BVH::T_tri / node.aabb.surfaceArea();
      for (uintp_t i = size - 1; i >= 1; --i) {
        uintp_t f = *(start + i);

        aabb.expandBy(V[F(0, f)]);
        aabb.expandBy(V[F(1, f)]);
        aabb.expandBy(V[F(2, f)]);

        real_t left_area = left_areas[i - 1];
        real_t right_area = aabb.surfaceArea();
        uintp_t prims_left = i;
        uintp_t prims_right = size - i;

        real_t sah_cost =
            2.0 * BVH::T_aabb +
            tri_factor * (prims_left * left_area + prims_right * right_area);
        if (sah_cost < best_cost) {
          best_cost = sah_cost;
          best_index = i;
          best_axis = axis;
        }
      }
    }

    if (best_index == -1) {
      /* Splitting does not reduce the cost, make a leaf */
      node.leaf.flag = 1;
      node.leaf.start = start - bvh.mIndices.data();
      node.leaf.size = size;
      return;
    }

    std::sort(start, end, [&](uintp_t f1, uintp_t f2) {
      return ((V[F(0, f1)](best_axis) + V[F(1, f1)](best_axis) +
               V[F(2, f1)](best_axis)) <
              (V[F(0, f2)](best_axis) + V[F(1, f2)](best_axis) +
               V[F(2, f2)](best_axis)));
    });

    uintp_t left_count = best_index;
    int node_idx_left = node_idx + 1;
    int node_idx_right = node_idx + 2 * left_count;
    node.inner.rightChild = node_idx_right;
    node.inner.unused = 0;

    execute_serially(bvh, node_idx_left, start, start + left_count, temp);
    execute_serially(bvh, node_idx_right, start + left_count, end,
                     temp + left_count);
  }
};

BVH::BVH(std::vector<uint32_t> *F, std::vector<Vector3f> *V,
         std::vector<Vector3f> *N, const AABB &aabb)
    : mF(F), mV(V), mN(N) {
  mNodes.resize(2 * mF->size() / 3);
  memset(mNodes.data(), 0, sizeof(BVHNode) * mNodes.size());
  mNodes[0].aabb = aabb;
  mIndices.resize(mF->size() / 3);
}

void BVH::build() {
  if (mF->size() == 0)
    return;

#if defined(SINGLE_PRECISION)
  if (sizeof(BVHNode) != 32)
    throw std::runtime_error(
        "BVH Node is not packed! Investigate compiler settings.");
#endif

  uintp_t total_size = mF->size() / 3;

  for (uintp_t i = 0; i < total_size; ++i)
    mIndices[i] = i;

  uintp_t *indices = mIndices.data();
#if defined(SINGLE_PRECISION)
  uintp_t *temp = new uintp_t[total_size];
  BVHBuildTask::execute_serially(*this, 0u, indices, indices + total_size, temp);
  delete[] temp;
#else
  real_t *temp = new real_t[total_size];
  BVHBuildTask::execute_serially(*this, 0u, indices, indices + total_size, temp);
  delete[] temp;
#endif  

  std::pair<real_t, uintp_t> stats = statistics();

  std::vector<BVHNode> compressed(stats.second);
  std::vector<uintp_t> skipped_accum(mNodes.size());

  for (int64_t i = stats.second - 1, j = mNodes.size(), skipped = 0; i >= 0;
       --i) {
    while (mNodes[--j].isUnused())
      skipped++;
    BVHNode &new_node = compressed[i];
    new_node = mNodes[j];
    skipped_accum[j] = skipped;

    if (new_node.isInner()) {
      new_node.inner.rightChild =
          i + new_node.inner.rightChild - j -
          (skipped - skipped_accum[new_node.inner.rightChild]);
    }
  }

  mNodes = std::move(compressed);
}

bool BVH::rayIntersect(Ray ray, uint32_t &idx, real_t &t, Vector2f *uv) const {
  if (mNodes.empty())
    return false;

  uintp_t node_idx = 0, stack[64];
  uintp_t stack_idx = 0;
  bool hit = false;
  t = std::numeric_limits<real_t>::infinity();

  while (true) {
    const BVHNode &node = mNodes[node_idx];

    if (!node.aabb.rayIntersect(ray)) {
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
      continue;
    }

    if (node.isInner()) {
      stack[stack_idx++] = node.inner.rightChild;
      node_idx++;
      assert(stack_idx < 64);
    } else {
      real_t _t;
      Vector2f _uv;
      for (uintp_t i = node.start(), end = node.end(); i < end; ++i) {
        if (rayIntersectTri(ray, mIndices[i], _t, _uv)) {
          idx = mIndices[i];
          t = ray.maxt = _t;
          hit = true;
          if (uv)
            *uv = _uv;
        }
      }
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
      continue;
    }
  }

  return hit;
}

bool BVH::rayIntersect(Ray ray, std::vector<uint32_t> &idx,
                       std::vector<Vector2f> *uv) const {
  if (mNodes.empty())
    return false;

  uintp_t node_idx = 0, stack[64];
  uintp_t stack_idx = 0;
  bool hit = false;

  while (true) {
    const BVHNode &node = mNodes[node_idx];

    if (!node.aabb.rayIntersect(ray)) {
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
      continue;
    }

    if (node.isInner()) {
      stack[stack_idx++] = node.inner.rightChild;
      node_idx++;
      assert(stack_idx < 64);
    } else {
      real_t _t;
      Vector2f _uv;
      for (uintp_t i = node.start(), end = node.end(); i < end; ++i) {
        if (rayIntersectTri(ray, mIndices[i], _t, _uv)) {
          idx.push_back(mIndices[i]);
          if (uv)
            uv->push_back(_uv);
          hit = true;
        }
      }
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
    }
  }

  return hit;
}

bool BVH::rayIntersect(Ray ray) const {
  if (mNodes.empty())
    return false;

  uintp_t node_idx = 0, stack[64];
  uintp_t stack_idx = 0;

  while (true) {
    const BVHNode &node = mNodes[node_idx];

    if (!node.aabb.rayIntersect(ray)) {
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
      continue;
    }

    if (node.isInner()) {
      stack[stack_idx++] = node.inner.rightChild;
      node_idx++;
      assert(stack_idx < 64);
    } else {
      real_t t;
      Vector2f uv;
      for (uintp_t i = node.start(), end = node.end(); i < end; ++i)
        if (rayIntersectTri(ray, mIndices[i], t, uv))
          return true;
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
      continue;
    }
  }

  return false;
}

void BVH::findNearestWithRadius(const Vector3f &p, real_t radius,
                                std::vector<uint32_t> &result,
                                bool includeSelf) const {
  result.clear();

  uintp_t node_idx = 0, stack[64];
  uintp_t stack_idx = 0;
  real_t radius2 = radius * radius;

  while (true) {
    const BVHNode &node = mNodes[node_idx];
    if (node.aabb.squaredDistanceTo(p) > radius2) {
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
      continue;
    }

    if (node.isInner()) {
      stack[stack_idx++] = node.inner.rightChild;
      node_idx++;
      assert(stack_idx < 64);
    } else {
      uintp_t start = node.leaf.start, end = start + node.leaf.size;
      for (uintp_t i = start; i < end; ++i) {
        uintp_t f = mIndices[i];
        Vector3f pointPos = Vector3f::Zero();
        for (int j = 0; j < 3; ++j)
          pointPos += mV->at((*mF)[j + f * 3]);
        pointPos *= 1.0 / 3.0;
        real_t pointDist2 = (pointPos - p).squaredNorm();
        if (pointDist2 < radius2 && (pointDist2 != 0 || includeSelf))
          result.push_back(f);
      }
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
      continue;
    }
  }
}

uint32_t BVH::findNearest(const Vector3f &p, real_t &radius,
                          bool includeSelf) const {
  uintp_t node_idx = 0, stack[64];
  uintp_t stack_idx = 0;
  real_t radius2 = radius * radius;
  uintp_t result = (uintp_t)-1;

  while (true) {
    const BVHNode &node = mNodes[node_idx];
    if (node.aabb.squaredDistanceTo(p) > radius2) {
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
      continue;
    }

    if (node.isInner()) {
      uintp_t left = node_idx + 1, right = node.inner.rightChild;
      real_t distLeft = mNodes[left].aabb.squaredDistanceTo(p);
      real_t distRight = mNodes[right].aabb.squaredDistanceTo(p);
      if (distLeft < distRight) {
        node_idx = left;
        if (distRight < radius2)
          stack[stack_idx++] = right;
      } else {
        node_idx = right;
        if (distLeft < radius2)
          stack[stack_idx++] = left;
      }
      assert(stack_idx < 64);
    } else {
      uintp_t start = node.leaf.start, end = start + node.leaf.size;
      for (uintp_t i = start; i < end; ++i) {
        uintp_t f = mIndices[i];
        Vector3f pointPos = Vector3f::Zero();
        for (int j = 0; j < 3; ++j)
          pointPos += mV->at((*mF)[j + f * 3]);
        pointPos *= 1.0 / 3.0;
        real_t pointDist2 = (pointPos - p).squaredNorm();

        if (pointDist2 < radius2 && (pointDist2 != 0 || includeSelf)) {
          radius2 = pointDist2;
          result = f;
        }
      }
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
      continue;
    }
  }
  radius = std::sqrt(radius2);
  return result;
}

void BVH::findKNearest(const Vector3f &p, uint32_t k, real_t &radius,
                       std::vector<std::pair<real_t, uint32_t>> &result,
                       bool includeSelf) const {
  result.clear();

  uintp_t node_idx = 0, stack[64];
  uintp_t stack_idx = 0;
  real_t radius2 = radius * radius;
  bool isHeap = false;
  auto comp = [](const std::pair<real_t, uintp_t> &v1,
                 const std::pair<real_t, uintp_t> &v2) {
    return v1.first < v2.first;
  };

  while (true) {
    const BVHNode &node = mNodes[node_idx];
    if (node.aabb.squaredDistanceTo(p) > radius2) {
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
      continue;
    }

    if (node.isInner()) {
      uintp_t left = node_idx + 1, right = node.inner.rightChild;
      real_t distLeft = mNodes[left].aabb.squaredDistanceTo(p);
      real_t distRight = mNodes[right].aabb.squaredDistanceTo(p);
      if (distLeft < distRight) {
        node_idx = left;
        if (distRight < radius2)
          stack[stack_idx++] = right;
      } else {
        node_idx = right;
        if (distLeft < radius2)
          stack[stack_idx++] = left;
      }
      assert(stack_idx < 64);
    } else {
      uintp_t start = node.leaf.start, end = start + node.leaf.size;
      for (uintp_t i = start; i < end; ++i) {
        uintp_t f = mIndices[i];
        Vector3f pointPos = Vector3f::Zero();
        for (int j = 0; j < 3; ++j)
          pointPos += mV->at((*mF)[j * 3 + f]);
        pointPos *= 1.0 / 3.0;
        real_t pointDist2 = (pointPos - p).squaredNorm();

        if (pointDist2 < radius2 && (pointDist2 != 0 || includeSelf)) {
          if (result.size() < k) {
            result.push_back(std::make_pair(pointDist2, f));
          } else {
            if (!isHeap) {
              /* Establish the max-heap property */
              std::make_heap(result.begin(), result.end(), comp);
              isHeap = true;
            }

            result.push_back(std::make_pair(pointDist2, f));
            std::push_heap(result.begin(), result.end(), comp);
            std::pop_heap(result.begin(), result.end(), comp);
            result.pop_back();

            /* Reduce the search radius accordingly */
            radius2 = result[0].first;
          }
        }
      }
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
      continue;
    }
  }
  radius = std::sqrt(radius2);
}

void BVH::findKNearest(const Vector3f &p, const Vector3f &n, uint32_t k,
                       real_t &radius,
                       std::vector<std::pair<real_t, uint32_t>> &result,
                       real_t angleThresh, bool includeSelf) const {
  result.clear();

  uintp_t node_idx = 0, stack[64];
  uintp_t stack_idx = 0;
  real_t radius2 = radius * radius;
  bool isHeap = false;
  angleThresh = std::cos(angleThresh * M_PI / 180);
  auto comp = [](const std::pair<real_t, uintp_t> &v1,
                 const std::pair<real_t, uintp_t> &v2) {
    return v1.first < v2.first;
  };

  while (true) {
    const BVHNode &node = mNodes[node_idx];
    if (node.aabb.squaredDistanceTo(p) > radius2) {
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
      continue;
    }

    if (node.isInner()) {
      uintp_t left = node_idx + 1, right = node.inner.rightChild;
      real_t distLeft = mNodes[left].aabb.squaredDistanceTo(p);
      real_t distRight = mNodes[right].aabb.squaredDistanceTo(p);
      if (distLeft < distRight) {
        node_idx = left;
        if (distRight < radius2)
          stack[stack_idx++] = right;
      } else {
        node_idx = right;
        if (distLeft < radius2)
          stack[stack_idx++] = left;
      }
      assert(stack_idx < 64);
    } else {
      uintp_t start = node.leaf.start, end = start + node.leaf.size;
      for (uintp_t i = start; i < end; ++i) {
        uintp_t f = mIndices[i];
        Vector3f pointPos = Vector3f::Zero();
        for (int j = 0; j < 3; ++j)
          pointPos += mV->at((*mF)[j + 3 * f]);
        pointPos *= 1.0 / 3.0;
        Vector3f pointNormal = Vector3f::Zero();
        for (int j = 0; j < 3; ++j)
          pointNormal += mN->at((*mF)[j + 3 * f]);
        real_t pointDist2 = (pointPos - p).squaredNorm();

        if (pointDist2 < radius2 && (pointDist2 != 0 || includeSelf) &&
            pointNormal.dot(n) > angleThresh) {
          if (result.size() < k) {
            result.push_back(std::make_pair(pointDist2, f));
          } else {
            if (!isHeap) {
              /* Establish the max-heap property */
              std::make_heap(result.begin(), result.end(), comp);
              isHeap = true;
            }

            result.push_back(std::make_pair(pointDist2, f));
            std::push_heap(result.begin(), result.end(), comp);
            std::pop_heap(result.begin(), result.end(), comp);
            result.pop_back();

            /* Reduce the search radius accordingly */
            radius2 = result[0].first;
          }
        }
      }
      if (stack_idx == 0)
        break;
      node_idx = stack[--stack_idx];
      continue;
    }
  }
  radius = std::sqrt(radius2);
}

bool BVH::rayIntersectTri(const Ray &ray, uint32_t i, real_t &t,
                          Vector2f &uv) const {
  const Vector3f &p0 = mV->at((*mF)[i * 3]), &p1 = mV->at((*mF)[1 + 3 * i]),
                 &p2 = mV->at((*mF)[2 + 3 * i]);

  Vector3f edge1 = p1 - p0, edge2 = p2 - p0;
  Vector3f pvec = ray.d.cross(edge2);

  real_t det = edge1.dot(pvec);
  if (det == 0.0)
    return false;
  real_t inv_det = 1.0 / det;

  Vector3f tvec = ray.o - p0;
  real_t u = tvec.dot(pvec) * inv_det;
  if (u < 0.0 || u > 1.0)
    return false;

  Vector3f qvec = tvec.cross(edge1);
  real_t v = ray.d.dot(qvec) * inv_det;

  if (v < 0.0 || u + v > 1.0)
    return false;

  real_t tempT = edge2.dot(qvec) * inv_det;
  if (tempT < ray.mint || tempT > ray.maxt)
    return false;

  t = tempT;
  uv << u, v;
  return true;
}

void BVH::printStatistics() const {
  cout << endl;
  cout << "Bounding Volume Hierarchy statistics:" << endl;
  cout << "    Tree nodes         : "
       << memString(sizeof(BVHNode) * mNodes.size()) << endl;
  cout << "    Index buffer       : " << memString(sizeof(uintp_t) * mF->size())
       << endl;
  cout << "    Total              : "
       << memString(sizeof(BVHNode) * mNodes.size() +
                    sizeof(uintp_t) * mF->size())
       << endl;
}

std::pair<real_t, uintp_t> BVH::statistics(uintp_t node_idx) const {
  const BVHNode &node = mNodes[node_idx];
  if (node.isLeaf()) {
    return std::make_pair(T_tri * node.leaf.size, 1u);
  } else {
    std::pair<real_t, uintp_t> stats_left = statistics(node_idx + 1u);
    std::pair<real_t, uintp_t> stats_right = statistics(node.inner.rightChild);
    real_t saLeft = mNodes[node_idx + 1u].aabb.surfaceArea();
    real_t saRight = mNodes[node.inner.rightChild].aabb.surfaceArea();
    real_t saCur = node.aabb.surfaceArea();
    real_t sahCost =
        2 * BVH::T_aabb +
        (saLeft * stats_left.first + saRight * stats_right.first) / saCur;

    return std::make_pair(sahCost, stats_left.second + stats_right.second + 1u);
  }
}

BVH::~BVH() {}
