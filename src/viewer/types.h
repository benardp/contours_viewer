//=============================================================================
// Copyright (C) 2013 Graphics & Geometry Group, Bielefeld University
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Library General Public License
// as published by the Free Software Foundation, version 2.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//=============================================================================


#ifndef SURFACE_MESH_TYPES_H
#define SURFACE_MESH_TYPES_H


//== INCLUDES =================================================================


#include <common.h>


//=============================================================================
namespace Eigen {
  template<typename A,typename B>
  typename internal::traits<A>::Scalar
  dot(const MatrixBase<A>& a,const MatrixBase<B>& b) {
    return a.dot(b);
  }

  template<typename A,typename B>
  inline typename MatrixBase<A>::template cross_product_return_type<B>::type
  cross(const MatrixBase<A>& a,const MatrixBase<B>& b) {
    return a.cross(b);
  }
}

namespace surface_mesh {


//=============================================================================

typedef Vector3f Vec3f;
typedef Vector2f Vec2f;

/// Scalar type
typedef real_t Scalar;

/// Point type
typedef Vector3f Point;

/// 3D vector type
typedef Vector3f Vec3;

/// Normal type
typedef Vector3f Normal;

/// Color type
typedef Vector3f Color;

/// Texture coordinate type
typedef Vector3f Texture_coordinate;


//=============================================================================
} // namespace surface_mesh
//=============================================================================
#endif // SURFACE_MESH_TYPES_H
//============================================================================
