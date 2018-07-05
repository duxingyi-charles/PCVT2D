/*
 *  _____ _____ ________  _
 * /  __//  __//  __/\  \//
 * | |  _|  \  |  \   \  / 
 * | |_//|  /_ |  /_  /  \ 
 * \____\\____\\____\/__/\\
 *
 * Graphics Environment for EXperimentations.
 *  Copyright (C) 2006 INRIA - Project ALICE
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: 
 *
 *     ALICE Project - INRIA
 *     INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 *  Note that the GNU General Public License does not permit incorporating
 *  the Software into proprietary programs. 
 */


#ifndef __LLOYD_ENERGY__
#define __LLOYD_ENERGY__

#include "geometry.h"

namespace Geex {

    // Note: we do not divide by 30 (since we do not care),
    // but we should normally ...

    template <class T> inline T Lloyd_energy(
        const vec2g<T>& p1, const vec2g<T>& p2, const vec2g<T>& p3
    ) {
        vec2g<T> U = p2 - p1 ;
        vec2g<T> V = p3 - p1 ;
		///dxy change : divide 30
		return det(U,V) * (U.length2() + V.length2() + dot(U,V)) /30.0 ;
//         return det(U,V) * (
//             U.length2() + V.length2() + dot(U,V) 
//         ) ;
    }

	//dxy copied from CVTRemesh

	// weighted Lloyd energy of a triangle which belong to Voronoi region of seed
	// weight is defined at each vertex and linear interpolated within triangle
	template <class T> inline T Lloyd_energy(
		//const vec3g<T>& seed, //seed = p1
		const vec2g<T>& p1, const vec2g<T>& p2, const vec2g<T>& p3,
		T a, T b, T c,
		vec2g<T>& gradient/*, T& V*/
		) {

			T t_area = triangle_area(p1, p2, p3);
			vec2g<T> sp[3];
// 			sp[0] = seed - p1;
// 			sp[1] = seed - p2;
// 			sp[2] = seed - p3;
			//
			sp[0] = p1 - p1;
			sp[1] = p1 - p2;
			sp[2] = p1 - p3;
			//

			T rho[3];
			T alpha[3];
			T Sp = a + b + c;
			//V  = t_area * Sp / 3.;
			rho[0] = a; rho[1] = b; rho[2] = c;
			for (int i =0; i < 3; i++)
			{
				alpha[i] = Sp + rho[i];
			}
			T f = 0;
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j <= i; j++)
				{
					f += (alpha[i] + rho[j]) * dot(sp[i], sp[j]);
				}
			}
			gradient = (t_area/6.0) * (4.0*Sp*p1 - (alpha[0]*p1+alpha[1]*p2+alpha[2]*p3));
			return t_area * f / 30.0;
	}

	//dxy end
}


#endif
