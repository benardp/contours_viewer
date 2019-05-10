#include <cmath>

int intersect2dSeg2dSegParametric(double p1[2], double p2[2], double p3[2],
                                  double p4[2], double &t, double &u,
                                  double epsilon);

// Compute the intersection between the 2D segments (p1,p2) and (p3,p4)
// returns 0 if they don't intersection,
//         1 if they do intersect: t and u then contain the parametric
//         locations of the intersection point along the two segments, 
//         2 if they are colinear.
int intersect2dSeg2dSegParametric(double p1[2], double p2[2], double p3[2],
                                  double p4[2], double &t, double &u,
                                  double epsilon) {
  double a1, a2, b1, b2, c1, c2; // Coefficients of line eqns
  double r1, r2, r3, r4;         // 'Sign' values
  double denom, num;             // Intermediate values

  // Compute a1, b1, c1, where line joining points p1 and p2
  // is "a1 x  +  b1 y  +  c1  =  0".
  a1 = p2[1] - p1[1];
  b1 = p1[0] - p2[0];
  c1 = p2[0] * p1[1] - p1[0] * p2[1];

  // Compute r3 and r4.
  r3 = a1 * p3[0] + b1 * p3[1] + c1;
  r4 = a1 * p4[0] + b1 * p4[1] + c1;

  // Check signs of r3 and r4.  If both point 3 and point 4 lie on
  // same side of line 1, the line segments do not intersect.
  if (r3 != 0 && r4 != 0 && r3 * r4 > 0.0)
    return 0;

  // Compute a2, b2, c2
  a2 = p4[1] - p3[1];
  b2 = p3[0] - p4[0];
  c2 = p4[0] * p3[1] - p3[0] * p4[1];

  // Compute r1 and r2
  r1 = a2 * p1[0] + b2 * p1[1] + c2;
  r2 = a2 * p2[0] + b2 * p2[1] + c2;

  // Check signs of r1 and r2.  If both point 1 and point 2 lie
  // on same side of second line segment, the line segments do
  // not intersect.
  if (r1 != 0 && r2 != 0 && r1 * r2 > 0.0)
    return 0;

  // Line segments intersect: compute intersection point.
  denom = a1 * b2 - a2 * b1;
  if (std::abs(denom) < epsilon)
    return 2;

  double d1, d2, e1;

  d1 = p1[1] - p3[1];
  d2 = p2[1] - p1[1];
  e1 = p1[0] - p3[0];

  num = -b2 * d1 - a2 * e1;
  t = num / denom;

  num = -b1 * d1 - a1 * e1;
  u = num / denom;

  return 1;
}
