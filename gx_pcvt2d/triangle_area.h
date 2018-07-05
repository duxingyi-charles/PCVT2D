// triangle_area.h - 

// Example code for "Efficient Generatino of Poisson-Disk Sampling
// Patterns," Thouis R. Jones, JGT vol. 11, No. 2, pp. 27-36
// 
// Copyright 2004-2006, Thouis R. Jones
// This code is distributed under the terms of the LGPL.


#ifndef __TRIANGLE_AREA_H__
#define __TRIANGLE_AREA_H__

#include "geometry.h"

namespace Geex {

double triangle_area_exclude(const vec2& p0, const vec2& p1, const vec2& p2, double exclude_radius);
vec2   random_point_in_triangle_exclude(const vec2& v1, const vec2& v2, const vec2& v3, double exclude_radius);

} // end of namespace Geex

//double triangle_area_exclude(GtsTriangle *t, GtsVertex *v, double exclude_radius);
//GtsVertex *random_point_in_triangle_exclude(GtsTriangle *t, GtsVertex *v, double exclude_radius);

#endif 