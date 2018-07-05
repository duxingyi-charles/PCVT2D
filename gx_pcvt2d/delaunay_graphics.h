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


#ifndef __DELAUNAY_GRAPHICS__
#define __DELAUNAY_GRAPHICS__

#include "delaunay.h"
#include <Geex/graphics/opengl.h>

namespace Geex {

    class DelaunayCVT ;
	class DelaunayIO ;
    
    class DelaunayGraphics {
    public:
        DelaunayGraphics(DelaunayCVT* CVT, DelaunayIO* IO=nil) ;
        void draw() ;

        GLboolean& show_domain() { return show_domain_ ; }
		GLboolean& show_vertices() { return show_vertices_; }
        GLfloat& vertices_size() { return vertices_size_ ; }
		GLboolean& vertices_color() { return vertices_color_ ; }
        GLfloat& centers_size() { return centers_size_ ; }
        GLboolean& show_primal_mesh() { return show_primal_mesh_ ; }
        GLboolean& show_dual_mesh() { return show_dual_mesh_ ; }
        GLboolean& colorize() { return colorize_ ; }
        GLboolean& show_cells() { return show_cells_ ; }
        GLboolean& show_energy() { return show_energy_ ; }
        GLboolean& show_field() { return show_field_ ; }
        GLfloat& quad_ratio() { return quad_ratio_ ; }
		GLboolean& show_boundary_cells() { return show_boundary_cells_ ; }
		GLboolean& show_pvd_euclidean() { return show_pvd_euclidean_ ; }
		GLboolean& show_copies() { return show_copies_ ; }
		GLboolean& show_disk() { return show_disk_ ; }
		GLboolean& show_mesh() { return show_mesh_ ; }
		GLboolean& show_min_max() { return show_min_max_ ;  } 
		GLboolean& show_inner_voronoi() { return show_inner_voronoi_ ; } 
		GLboolean& show_regularity() { return show_regularity_ ; }

		///dxy add: getter
		GLboolean& show_regular_grad_only() { return show_regular_grad_only_; }
		GLboolean& show_lloyd_grad() { return show_lloyd_grad_; }
		GLboolean& show_direction_grad() { return show_direction_grad_ ; }
		GLboolean& show_total_grad() { return show_total_grad_; }
		GLboolean& show_constrained_edge() { return show_constrained_edge_; }
		GLboolean& show_direction_field() { return show_direction_field_; }
		GLboolean& show_relative_area() { return show_relative_area_; }
		GLboolean& show_triangle_area() { return show_triangle_area_; }
		GLboolean& show_long_edge() { return show_long_edge_; }
		GLboolean& show_short_edge() { return show_short_edge_; }
		GLboolean& show_edge_stretch() { return show_edge_stretch_; }
		GLboolean& show_separatrix() { return show_separatrix_ ;}
		GLboolean& show_selected() { return show_selected_; }

		int& lloyd_grad_magnification() { return lloyd_grad_magnification_; }
		int& direct_grad_magnification() { return direct_grad_magnification_; }
		int& total_grad_magnification() { return total_grad_magnification_; }
		int& stretch_long_rank() { return stretch_long_rank_; }
		int& stretch_short_rank() { return stretch_short_rank_; }
		int& stretch_all_rank() { return stretch_all_rank_; }
		///

    protected:
        void draw_domain() ;
        void draw_vertices() ;
        void draw_centers() ;
        void draw_cells(bool colorize = true, bool mesh = false) ;
        void draw_primal_mesh() ;
        void draw_dual_mesh() ;
		void draw_dual_edges() ;
        void draw_non_hex_cells() ;
		void draw_non_hex_periodic_cells() ;
        void draw_energy() ;
        void draw_field() ;
		void draw_boundary_cells() ;
		void draw_edge_hist() ;
		void draw_disk() ;
		void draw_mesh() ;
		void draw_selection() ;
		void draw_min_max() ;
		void draw_inner_voronoi() ;
        void is_min_max(Delaunay::Vertex_handle v, bool& is_min, bool& is_max) ;
        void draw_dual_facet(Delaunay::Vertex_handle v) ;
        double radius(Delaunay::Vertex_handle v) ;
		void draw_regularity() ;

		///dxy change
		/*void draw_non_hex_cells() ;
		void draw_energy() ;*/
// 		void draw_non_hex_cells(bool info = false) ;
// 		void draw_energy() ;
		///dxy add
		void draw_relative_area();
		void draw_primal_area();
		//
		void draw_lloyd_grad();
		void draw_direction_grad();
		void draw_regular_direction_grad();
		void draw_total_grad();
		void draw_direction_field();
		//
		void draw_constrained_edge();
		void draw_long_edge(); //7-7
		void draw_short_edge(); //5-5
		void draw_isolate_long_edge();
		void draw_isolate_short_edge();
		void draw_edge_stretch();
		//
		void draw_edge_stretch_period();
		//
		void draw_separatrix();
		void draw_separatrix_period();
		void draw_selected_separatrix();
		//
		void draw_selected();
		//auxiliary
		void rotate_vector(vec2& v, double angle);
		void draw_vector(const vec2&, const vec2&);
		void draw_primal_triangle(Delaunay::Face_handle f);
		void draw_triEdge(const triEdge&);
		void draw_dirEdge(const dirEdge&);
		void draw_Separatrix(const Separatrix&);
		void draw_Separatrix(const Separatrix&, float);
		//pcvt add
		void draw_period_dual_facet(Delaunay::Vertex_handle v);
		///

    private:
        Delaunay* delaunay_ ;
		DelaunayIO* IO_ ;
        DelaunayCVT* CVT_ ;
        GLboolean show_domain_ ;
		GLboolean show_vertices_;
        GLfloat vertices_size_ ;
		GLboolean vertices_color_ ;
        GLfloat centers_size_ ;
        GLboolean show_primal_mesh_ ;
        GLboolean show_dual_mesh_ ;
        GLboolean colorize_ ;
        GLboolean show_cells_ ;
        GLboolean show_energy_ ;
        GLboolean show_field_ ;
        bool non_convex_ ;
        GLfloat quad_ratio_ ;
		GLboolean show_boundary_cells_ ;
		GLboolean show_pvd_euclidean_ ;
		GLboolean show_copies_ ;
		//GLboolean show_edge_hist_ ;
		GLboolean show_disk_ ;
		GLboolean show_mesh_ ;
		GLboolean show_min_max_ ;
		GLboolean show_inner_voronoi_ ;
		GLboolean show_regularity_ ; 

		///dxy add
		GLboolean show_regular_grad_only_;
		GLboolean show_lloyd_grad_ ;
		GLboolean show_direction_grad_ ;
		GLboolean show_total_grad_;
		GLboolean show_constrained_edge_ ;
		GLboolean show_direction_field_;
		GLboolean show_relative_area_;
		GLboolean show_triangle_area_;
		GLboolean show_long_edge_;
		GLboolean show_short_edge_;
		GLboolean show_edge_stretch_;
		GLboolean show_separatrix_;
		GLboolean show_selected_;
		int lloyd_grad_magnification_;
		int direct_grad_magnification_;
		int total_grad_magnification_;
		int stretch_long_rank_;
		int stretch_short_rank_;
		int stretch_all_rank_;
		///
	} ;

} 

#endif
