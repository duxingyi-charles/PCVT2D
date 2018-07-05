/* triangle_area.c - */

// Example code for "Efficient Generatino of Poisson-Disk Sampling
// Patterns," Thouis R. Jones, JGT vol. 11, No. 2, pp. 27-36
// 
// Copyright 2004-2006, Thouis R. Jones
// This code is distributed under the terms of the LGPL.

#include "triangle_area.h"

//#include <stdlib.h>
//#include <cmath>
//#include "gts.h"
#define printf
//#define g_assert(x) 

namespace Geex {

	double sgn(double a)
	{
	  if (a < 0.0) return (-1.0);
	  return (1.0);
	}

	int between(double a, double b, double c)
	{
	return (((a <= b) && (b <= c)) ||
			((a >= b) && (b >= c)));
	}

	void clip_to_circle_origin(double r, double x1, double y1, double x2, double y2,
							double *ix1, double *iy1, double *ix2, double *iy2)
	{
		double dx, dy, dr2, D;

		// from MathWorld

		dx = x2 - x1;
		dy = y2 - y1;
		dr2 = dx*dx + dy*dy;
		D = x1 * y2 - x2 * y1;

		if ((r*r*dr2-D*D) <= 0.0) {
		// no intersection?
			std::cerr << "WARNING: no intersection in clip_to_circle_origin\n";
			*ix1 = x2; *iy1 = y2;
			return;
		}

		*ix1 = (D*dy + sgn(dy)*dx*sqrt(r*r*dr2-D*D))/dr2;
		*iy1 = (-D*dx + fabs(dy)*sqrt(r*r*dr2-D*D))/dr2;
  
		*ix2 = (D*dy - sgn(dy)*dx*sqrt(r*r*dr2-D*D))/dr2;
		*iy2 = (-D*dx - fabs(dy)*sqrt(r*r*dr2-D*D))/dr2;

		return;
	}

	#define SQR(a) ((a)*(a))

	vec2 clip_to_circle(const vec2& center, double exclude_radius, const vec2& out, const vec2& in)
	{
	  double x1, y1, x2, y2;
	  clip_to_circle_origin(exclude_radius, 
							out.x - center.x, out.y - center.y,
							in.x - center.x, in.y - center.y,
							&x1, &y1, &x2, &y2);

	  x1 += center.x;
	  y1 += center.y;
	  x2 += center.x;
	  y2 += center.y;

	  //printf("clipping %f %f     %f %f   cetner %f %f\n     %f %f      %f %f\n",
			// GTS_POINT(out)->x, GTS_POINT(out)->y, 
			// GTS_POINT(in)->x, GTS_POINT(in)->y, 
			// GTS_POINT(center)->x, GTS_POINT(center)->y,
			// x1, y1, x2, y2);
         

	  if ((SQR(x1 - out.x) + SQR(y1 - out.y)) <
		  (SQR(x2 - out.x) + SQR(y2 - out.y))) {

		return (vec2(x1, y1));
	  } else {
		return (vec2(x2, y2));
	  }
	}

	void closest_approach_circle_origin(double x1, double y1, double x2, double y2,
										double *ix, double *iy)
	{
	  double dx, dy, dr2, D;

	  // same as intersect above, but assume that discriminant is zero

	  dx = x2 - x1;
	  dy = y2 - y1;
	  dr2 = dx*dx + dy*dy;
	  D = x1 * y2 - x2 * y1;

	  *ix = (D*dy)/dr2;
	  *iy = (-D*dx)/dr2;

	  if (! ((between(x1, *ix, x2) && (between(y1, *iy, y2))))) {
		if ((SQR(x1-*ix)+SQR(y1-*iy)) < (SQR(x2-*ix)+SQR(y2-*iy))) {
		  *ix = x1;
		  *iy = y1;
		} else {
		  *ix = x2;
		  *iy = y2;
		}
	  }

	  return;
	}

	vec2 closest_approach_circle(const vec2& center, const vec2& v1, const vec2& v2)
	{
	  double x, y;
	  closest_approach_circle_origin(v1.x - center.x, v1.y - center.y,
									 v2.x - center.x, v2.y - center.y,
									 &x, &y);

	  return vec2(x + center.x, y + center.y);
	}

	void clip_to_circle_both(const vec2& center, double exclude_radius, const vec2& v1, const vec2& v2, vec2& out1, vec2& out2)
	{
	  vec2 closest_approach = closest_approach_circle(center, v1, v2);

	  gx_assert(distance(center, closest_approach) < exclude_radius);

	  out1 = clip_to_circle(center, exclude_radius, v1, closest_approach);
	  out2 = clip_to_circle(center, exclude_radius, v2, closest_approach);
   }

	int intersects_circle(const vec2& center, double exclude_radius, const vec2& v1, const vec2& v2)
	{
	  const vec2& closest_approach = closest_approach_circle(center, v1, v2);

	  if (distance(center, closest_approach) < exclude_radius) {
		return (1);
	  } else {
		return (0);
	  }
	}

	double threepts_area_exclude_compute(const vec2& v1, const vec2& v2, const vec2& v3, double exclude_radius)
	{
	  vec2 v12, v13;
	  double triarea, angle, circarea, z;

	  v12 = v2 - v1 ; //gts_vector_init(v12, GTS_POINT(v1), GTS_POINT(v2));
	  v13 = v3 - v1 ; //gts_vector_init(v13, GTS_POINT(v1), GTS_POINT(v3));

	  z = fabs(v12[0] * v13[1] - v12[1] * v13[0]);
	  triarea = z / 2.0;
  
	  if (z == 0.0) {
		  angle = 0.0;
	  } else {
		  double valsin = z / (v12.length() * v13.length()) ;
		  if(valsin > 1.0) valsin = 1.0 ;
		  angle = asin(valsin) ; //v12.length() * v13.length()));
	  }

	  circarea = 0.5 * angle * exclude_radius * exclude_radius;

	  if ((triarea - circarea) <= 0.0) {
		return (0.0);
	  } else {
		return (triarea - circarea);
	  }
	}

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! v1 mast be the fisrt point !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	double threepts_area_exclude(const vec2& v1, const vec2& v2, const vec2& v3, double exclude_radius)
	{
	  int intersects12 = (distance(v1, v2) > exclude_radius);
	  int intersects13 = (distance(v1, v3) > exclude_radius);


	  // easy case
	  if ((! intersects12) & (! intersects13)) return (0.0);

	  // one vertex out
	  if ((intersects12) & (! intersects13)) {
		vec2 v23 = clip_to_circle(v1, exclude_radius, v2, v3);
		double area = threepts_area_exclude_compute(v1, v2, v23, exclude_radius);
		return area;
	  }

	  // one vertex out
	  if ((! intersects12) & (intersects13)) {
		vec2 v32 = clip_to_circle(v1, exclude_radius, v3, v2);
		double area = threepts_area_exclude_compute(v1, v32, v3, exclude_radius);
		return area;
	  }

	  // two vertices out
	  if (intersects_circle(v1, exclude_radius, v2, v3)) {
		// two separate areas
		vec2 v23_near2, v32_near3;
		double area21, area31;

		clip_to_circle_both(v1, exclude_radius, v2, v3, v23_near2, v32_near3);

		area21 = threepts_area_exclude_compute(v1, v2, v23_near2, exclude_radius);
		area31 = threepts_area_exclude_compute(v1, v3, v32_near3, exclude_radius);
    
		return (area21 + area31);
	  } else {
		// one area
		return (threepts_area_exclude_compute(v1, v2, v3, exclude_radius));
	  }
	}

	double triangle_area_exclude(const vec2& p0, const vec2& p1, const vec2& p2, double exclude_radius) {
		return threepts_area_exclude(p0, p1, p2, exclude_radius) ;
	}

	static int same_point(const vec2& v1, const vec2& v2)
	{
	  return ((v1.x == v2.x) &&
			  (v1.y == v2.y));
	}

	void random_point_in_threepts(const vec2& v1,const vec2& v2, const vec2& v3, vec2& out)
	{
	  double a1 = Numeric::random_float64(), a2 = Numeric::random_float64();
  
	  while ((a1 + a2) > 1.0) {
		a1 = Numeric::random_float64(); a2 = Numeric::random_float64();
	  }
  
	  printf("a1, a2, %f %f\n", a1, a2);

	  printf("source: %f %f    %f %f    %f %f\n",
			 v1.x, v1.y, 
			 v2.x, v2.y, 
			 v3.x, v3.y);

	  out = vec2(	a1 * v1.x + a2 * v2.x + (1 - a1 - a2) * v3.x,
					a1 * v1.y + a2 * v2.y + (1 - a1 - a2) * v3.y ) ;
	}

	#define MAXCOUNT 100000
	vec2 random_point_in_threepts_exclude_compute(const vec2& center, const vec2& v1, const vec2& v2, const vec2& v3, double exclude_radius)
	{
	  int count = 0;
	  vec2 out(0.0, 0.0);


	  printf("dists: %e %e %e\n",
			 distance(v1, center) - exclude_radius,
			 distance(v2, center) - exclude_radius,
			 distance(v3, center) - exclude_radius);

	  gx_assert(distance(v1, center) >= (exclude_radius - 1e-15));
	  gx_assert(distance(v2, center) >= (exclude_radius - 1e-15));
	  gx_assert(distance(v3, center) >= (exclude_radius - 1e-15));



	  for (count = 0; count < MAXCOUNT; count++) {
	  printf("source: %f %f    %f %f    %f %f    (exc %f %f, %f)\n",
			 v1.x, v1.y, 
			 v2.x, v2.y, 
			 v3.x, v3.y,
			 center.x, center.y, exclude_radius);

	  printf("    Dists: %f %f %f\n",
			 distance(v1, center),
			 distance(v2, center),
			 distance(v3, center));



		random_point_in_threepts(v1, v2, v3, out);

		if (distance(center, out) >= exclude_radius)
		  return (out);
	  }

	  //fprintf(stderr, "Died\n");
	  std::cerr << "Died\n" ;

	  exit(-1);

	  return (out);
	}



	vec2 random_point_in_threepts_exclude(const vec2& v1, const vec2& v2, const vec2& v3, double exclude_radius)
	{
	  int intersects12 = (distance(v1, v2) > exclude_radius);
	  int intersects13 = (distance(v1, v3) > exclude_radius);

	  // one edge of triangle must be outside
	  gx_assert(intersects12 || intersects13);

	  // one vertex out -> v2
	  if ((intersects12) & (! intersects13)) {
		vec2 v21 = clip_to_circle(v1, exclude_radius, v2, v1);
		vec2 v23 = clip_to_circle(v1, exclude_radius, v2, v3);
		vec2 out = random_point_in_threepts_exclude_compute(v1, v2, v21, v23, exclude_radius);
		printf ("v2 out\n");
		return out;
	  }

	  // one vertex out -> v3
	  if ((! intersects12) & (intersects13)) {
		vec2 v31 = clip_to_circle(v1, exclude_radius, v3, v1);
		vec2 v32 = clip_to_circle(v1, exclude_radius, v3, v2);
		vec2 out = random_point_in_threepts_exclude_compute(v1, v3, v31, v32, exclude_radius);
		printf ("v3 out\n");
		return out;
	  }

	  // two vertices out 
	  if (intersects_circle(v1, exclude_radius, v2, v3)) {
		// two separate areas
		vec2 out;
		vec2 v23_near2, v32_near3;
		double area21, area31;

		printf ("both out - intersecting\n");

		clip_to_circle_both(v1, exclude_radius, v2, v3, v23_near2, v32_near3);

		area21 = threepts_area_exclude_compute(v1, v2, v23_near2, exclude_radius);
		area31 = threepts_area_exclude_compute(v1, v3, v32_near3, exclude_radius);
		if ((Numeric::random_float64() * (area21 + area31)) > area31) {
		  vec2 v21 = clip_to_circle(v1, exclude_radius, v2, v1);
		  out = random_point_in_threepts_exclude_compute(v1, v21, v23_near2, v2, exclude_radius);
		} else {
		  vec2 v31 = clip_to_circle(v1, exclude_radius, v3, v1);
		  out = random_point_in_threepts_exclude_compute(v1, v31, v32_near3, v3, exclude_radius);
		}
		return (out);
	  }

	  // final case, v2->v3 doesn't intersect circle
	  {
		vec2 out;
		vec2 closest_approach = closest_approach_circle(v1, v2, v3);
		double area12c = threepts_area_exclude_compute(v1, v2, closest_approach, exclude_radius);
		double area13c = threepts_area_exclude_compute(v1, v3, closest_approach, exclude_radius);
		printf ("both out - non-intersecting\n");

		if ((Numeric::random_float64() * (area12c + area13c)) > area13c && area13c!=0) {
		  // select from triangle 1,2,closest_approach (with exclusion)
		  vec2 vc1 = clip_to_circle(v1, exclude_radius, closest_approach, v1);
		  double area_all_out = fabs(triangle_area(v2, vc1, closest_approach));
		  if ((Numeric::random_float64() * area12c) > area_all_out) {
			vec2 v21 = clip_to_circle(v1, exclude_radius, v2, v1);
			// select from triangle v21, vc1, v2 (with exclusion)
			printf("case 1\n");
			printf("vc1 dist %f\n", distance(vc1, v1));
			printf("v2 dist %f\n", distance(v2, v1));
			printf("v21 dist %f\n", distance(v21, v1));

			printf("center %f %f\n", v1.x, v1.y);
			printf("vc1 %f %f\n", vc1.x, vc1.y);
			printf("v2 %f %f\n", v2.x, v2.y);
			printf("v21 %f %f\n", v21.x, v21.y);
			printf("closest_approach %f %f\n", closest_approach.x, closest_approach.y);

			out = random_point_in_threepts_exclude_compute(v1, v21, vc1, v2, exclude_radius);
			printf("out %f %f\n", out.x, out.y);

			gx_assert(distance(v1, out) >= exclude_radius);
		  } else {
			// select from triangle v2, vc1, closest_approach (no exclusion)
			out = vec2(0.0, 0.0);
			printf("case 2\n");
			printf("vc1 dist %f\n", distance(vc1, v1));
			printf("v2 dist %f\n", distance(v2, v1));
			printf("closest_approach dist %f (should be >= %f)\n", distance(closest_approach, v1), exclude_radius);
			printf("center %f %f\n", v1.x, v1.y);
			printf("vc1 %f %f\n", vc1.x, vc1.y);
			printf("v2 %f %f\n", v2.x, v2.y);
			printf("closest_approach %f %f\n", closest_approach.x, closest_approach.y);

			random_point_in_threepts(vc1, v2, closest_approach, out);
			printf("out %f %f\n", out.x, out.y);

			if(distance(v1, out) < exclude_radius) {
				std::cerr << "v1 " << v1 << " out " << out << " distance " << distance(v1, out) << std::endl ;
				std::cerr << "area 12 " << area12c << " area 13c " << area13c << std::endl ;
			}

			gx_assert(distance(v1, out) >= exclude_radius); 

		  }
		} else {
		  // select from triangle 1,3,closest_approach (with exclusion)
		  vec2 vc1 = clip_to_circle(v1, exclude_radius, closest_approach, v1);
		  double area_all_out = fabs(triangle_area(v3, vc1, closest_approach));
		  if ((Numeric::random_float64() * area13c) > area_all_out) {
			// select from triangle v31, vc1, v3 (with exclusion)
			printf("case 3\n");
			vec2 v31 = clip_to_circle(v1, exclude_radius, v3, v1);

			printf("center %f %f\n", v1.x, v1.y);
			printf("vc1 %f %f\n", vc1.x, vc1.y);
			printf("v3 %f %f\n", v3.x, v3.y);
			printf("v31 %f %f\n", v31.x, v31.y);
			printf("closest_approach %f %f\n", closest_approach.x, closest_approach.y);


			printf("vc1 dist %f %e\n", distance(vc1, v1), exclude_radius - distance(vc1, v1));
			printf("v31 dist %f\n", distance(v31, v1));
			printf("v3 dist %f\n", distance(v3, v1));
			printf("closest_approach dist %f (should be >= %f)\n", distance(closest_approach, v1), exclude_radius);

			gx_assert(distance(v3, v1) >= exclude_radius);

			out = random_point_in_threepts_exclude_compute(v1, v31, v3, vc1, exclude_radius);
			gx_assert(distance(v1, out) >= exclude_radius); 
		  } else {
			// select from triangle v3, vc1, closest_approach (no exclusion)
			out = vec2(0.0, 0.0);
			printf("case 4\n");
			printf("closest_approach %f %f\n", closest_approach.x, closest_approach.y);

			printf("vc1 %f %f\n", vc1.x, vc1.y);

			random_point_in_threepts(v3, vc1, closest_approach, out);
			printf("out %f %f\n", out.x, out.y);
			gx_assert(distance(v1, out) >= exclude_radius); 
		  }
		}
		return (out);
	  }
	}
	
	vec2   random_point_in_triangle_exclude(const vec2& v1, const vec2& v2, const vec2& v3, double exclude_radius) {
		return random_point_in_threepts_exclude(v1, v2, v3, exclude_radius) ;
	}
	//vec2 random_point_in_triangle_exclude(GtsTriangle *t, GtsVertex *v, double exclude_radius)
	//{
	//  GtsVertex *v1, *v2, *v3, *vout;
 // 
	//  gts_triangle_vertices(t, &v1, &v2, &v3);
	//  if (same_point(v, v2)) {
	//	v2 = v1;
	//	v1 = v;
	//  }
	//  else if (same_point(v, v3)) {
	//	v3 = v1;
	//	v1 = v;
	//  }

	//  printf("TAKING from %f %f \n   %f %f \n    %f %f\n",
	//		 GTS_POINT(v1)->x, GTS_POINT(v1)->y, 
	//		 GTS_POINT(v2)->x, GTS_POINT(v2)->y, 
	//		 GTS_POINT(v3)->x, GTS_POINT(v3)->y);

	//  g_assert(same_point(v, v1));
 //   
	//  g_assert((distance(GTS_POINT(v1), GTS_POINT(v2)) >= exclude_radius) ||
	//		   (distance(GTS_POINT(v1), GTS_POINT(v3)) >= exclude_radius));

	//  vout = random_point_in_threepts_exclude(v1, v2, v3, exclude_radius);

	//  printf("dist: %f\n", distance(GTS_POINT(v1), GTS_POINT(vout)));
	//  g_assert(distance(GTS_POINT(v1), GTS_POINT(vout)) >= exclude_radius);

	//  return (vout);
	//}

} // end of namespace Geex



