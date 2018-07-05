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

#include "delaunay_graphics.h"
#include <delaunay_cvt.h>
#include <delaunay_io.h>
#include <Geex/combinatorics/map.h>
#include <glut_viewer/glut_viewer.h>

namespace Geex {


	static const double c1 = 0.35 ;
	static const double c2 = 0.5 ;
	static const double c3 = 1.0 ;

	static double color_table[12][3] = 
	{
		{c3, c2, c2},
		{c2, c3, c2},
		{c2, c2, c3},
		{c2, c3, c3},
		{c3, c2, c3},
		{c3, c3, c2},

		{c1, c2, c2},
		{c2, c1, c2},
		{c2, c2, c1},
		{c2, c1, c1},
		{c1, c2, c1},
		{c1, c1, c2}

	} ;

	static int random_color_index_ = 0 ;


	static void gl_random_color() {
		glColor3f(
			color_table[random_color_index_][0], 
			color_table[random_color_index_][1], 
			color_table[random_color_index_][2]
		) ;
		random_color_index_ = (random_color_index_ + 1) % 12 ;
	}

	static void gl_random_color(int index) {
		random_color_index_ = index % 12 ;
		gl_random_color() ;
	}

	static void gl_randomize_colors() {
		random_color_index_ = 0 ;
	}

	static void gl_vertex_color(int degree) {
		if(degree==6) {
			glColor3f(0, 1, 0) ;
		}
		else if(degree==7) {
			glColor3f(1.0, 0.5, 0.25) ;
		}
		else if(degree==5) {
			glColor3f(0.42, 0.55, 0.9) ;
		}
		else if(degree>7) {
			glColor3f(0.5, 0, 0) ;
		}
		else {
			glColor3f(0, 0, 0.5) ;
		}
	}

	///dxy add
	//hsv ([0,360], [0,1], [0,1]),  rgb ([0,1]...)
	static void HSVtoRGB(real& fR, real& fG, real& fB, real& fH, real& fS, real& fV) {
		real fC = fV * fS; // Chroma
		real fHPrime = fmod(fH / 60.0, 6);
		real fX = fC * (1 - fabs(fmod(fHPrime, 2) - 1));
		real fM = fV - fC;

		if(0 <= fHPrime && fHPrime < 1) {
			fR = fC;
			fG = fX;
			fB = 0;
		} else if(1 <= fHPrime && fHPrime < 2) {
			fR = fX;
			fG = fC;
			fB = 0;
		} else if(2 <= fHPrime && fHPrime < 3) {
			fR = 0;
			fG = fC;
			fB = fX;
		} else if(3 <= fHPrime && fHPrime < 4) {
			fR = 0;
			fG = fX;
			fB = fC;
		} else if(4 <= fHPrime && fHPrime < 5) {
			fR = fX;
			fG = 0;
			fB = fC;
		} else if(5 <= fHPrime && fHPrime < 6) {
			fR = fC;
			fG = 0;
			fB = fX;
		} else {
			fR = 0;
			fG = 0;
			fB = 0;
		}

		fR += fM;
		fG += fM;
		fB += fM;
	}
	static void hsv2rgb(vec3& hsv, vec3& rgb) 
	{
		real& fR = rgb[0];
		real& fG = rgb[1];
		real& fB = rgb[2];
		real& fH = hsv[0];
		real& fS = hsv[1];
		real& fV = hsv[2];
		HSVtoRGB(fR, fG, fB, fH, fS, fV);
	}
	///dxy add end

	DelaunayGraphics::DelaunayGraphics(DelaunayCVT* CVT, DelaunayIO* IO) 
		: delaunay_(CVT->delaunay())
		, CVT_(CVT)
		, IO_(IO) {
			show_domain_ = GL_TRUE ; 
			show_vertices_ = GL_TRUE;
			vertices_size_ = 0.1 ;
			vertices_color_ = GL_TRUE ;
			centers_size_ = 0.0 ;
			show_primal_mesh_ = GL_FALSE ;
			show_dual_mesh_ = GL_FALSE/*GL_TRUE*/ ;
			colorize_ = GL_FALSE ;
			show_cells_ = GL_FALSE ;
			show_energy_ = GL_FALSE ;
			show_field_ = GL_FALSE ;
			quad_ratio_ = 0.0 ;
			show_boundary_cells_ = GL_FALSE ;
			show_pvd_euclidean_ = GL_FALSE ;
			show_copies_ = GL_FALSE ;
			show_disk_ = GL_FALSE ;
			show_mesh_ = GL_FALSE ;
			show_min_max_ = GL_FALSE ;
			show_inner_voronoi_ = GL_FALSE ;
			show_regularity_ = GL_FALSE ;

			///dxy add
			show_regular_grad_only_ = GL_FALSE;
			show_lloyd_grad_ = GL_FALSE;
			show_direction_grad_ = GL_FALSE;
			show_constrained_edge_ = GL_FALSE;
			show_relative_area_ = GL_FALSE;
			show_triangle_area_ = GL_FALSE;
			show_separatrix_ = GL_FALSE;
			show_edge_stretch_ = GL_FALSE;
			show_direction_field_ = GL_FALSE;
			show_lloyd_grad_ = GL_FALSE;
			show_direction_field_ = GL_FALSE;
			show_total_grad_ = GL_FALSE;
			show_selected_ = GL_TRUE /*GL_FALSE*/;
			show_long_edge_ = GL_FALSE;
			show_short_edge_ = GL_FALSE;

			lloyd_grad_magnification_ = 1;
			direct_grad_magnification_ = 1;
			total_grad_magnification_ = 1;
			stretch_long_rank_  = -1;
			stretch_short_rank_ = -1;
			stretch_all_rank_ = -1;
			///
	}

	void DelaunayGraphics::draw() {
		non_convex_ = delaunay_->non_convex_mode_ ;
		glDisable(GL_LIGHTING) ;
		glUseProgramObjectARB(0) ;

		gl_randomize_colors() ;

		if(show_domain_) {
			draw_domain() ;
		}

		if(vertices_size_ != 0.0 && show_vertices_) {
			draw_vertices() ;
		}

		if(centers_size_ != 0.0) {
			draw_centers() ;
		}

		/*if(show_primal_mesh_) {
		draw_primal_mesh() ;
		}*/

		//if(show_inner_voronoi_) {
		//	draw_inner_voronoi() ;
		//}

		//if(show_dual_mesh_) {
		//	//            draw_dual_mesh() ;
		//	draw_dual_edges() ;
		//}



		if(show_cells_) {
			if(delaunay_->period())
				draw_non_hex_periodic_cells() ;
			else
				draw_non_hex_cells() ;
		}

		if(show_energy_) {
			draw_energy() ;
		}

		///dxy add

		if (show_relative_area_)
		{
			draw_relative_area();
		}
		if (show_triangle_area_)
		{
			draw_primal_area();
		}
		///



		if(show_field_) {
			draw_field() ;
		}

		if(show_boundary_cells_) {
			draw_boundary_cells() ;
		}

		if(show_disk_) {
			draw_disk() ;
		}

		if(show_min_max_) {
			draw_min_max() ;
		}

		if(show_mesh_) {
			draw_mesh() ;
		}

		if(IO_->show_edge_hist()) {
			draw_edge_hist() ;
		}

		if(show_regularity_) {
			draw_regularity() ;
		}

		draw_selection() ;

		///dxy add

		glDisable(GL_TEXTURE_2D) ;
		glDisable(GL_TEXTURE_1D) ;
		glUseProgramObjectARB(0) ;

		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST) ;
		glHint(GL_POINT_SMOOTH_HINT, GL_NICEST) ;
		glEnable(GL_LINE_SMOOTH) ;
		glEnable(GL_POINT_SMOOTH) ;
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA) ;
		glEnable(GL_BLEND) ;
		//glDepthMask(GL_FALSE) ;

		/*if (show_constrained_edge_)
		{
		draw_constrained_edge();
		}*/

		if (show_selected_)
		{
			draw_selected();
			draw_selected_separatrix();
		}

		if (show_lloyd_grad_)
		{
			draw_lloyd_grad();
		}

		if (show_direction_grad_)
		{
			if (show_regular_grad_only_)
			{
				draw_regular_direction_grad();
			}
			else draw_direction_grad();
		}

		if (show_total_grad_)
		{
			draw_total_grad();
		}

		if (show_direction_field_)
		{
			draw_direction_field();
		}

		if (show_long_edge_)
		{
			draw_isolate_long_edge();
		}
		if (show_short_edge_)
		{
			draw_isolate_short_edge();
		}

		if (show_edge_stretch_)
		{
			if (delaunay_->period())
			{
				draw_edge_stretch_period();
			}
			else {
				draw_edge_stretch();
			}
		}

		if (show_separatrix_)
		{
			if (delaunay_->period())
			{
				draw_separatrix_period();
			}
			else
			{
				draw_separatrix();
			}
		}



		///dxy move
		if(show_primal_mesh_) {
			draw_primal_mesh() ;
		}

		if(show_inner_voronoi_) {
			draw_inner_voronoi() ;
		}

		if(show_dual_mesh_) {
			draw_dual_edges() ;
		}
		///

		glDisable(GL_LINE_SMOOTH) ;
		glDisable(GL_BLEND) ;
		//glDepthMask(GL_TRUE) ;
		///
	}

	void DelaunayGraphics::draw_domain() {
		glLineWidth(3) ;
		glColor3f(0.1, 0.1, 0.1) ;
		glBegin(GL_LINES) ;
		for(unsigned int i=0; i<delaunay_->boundary_.size(); i++) {
			glVertex(delaunay_->boundary_[i].vertex[0]) ;
			glVertex(delaunay_->boundary_[i].vertex[1]) ;
		}
		glEnd() ;
	}

	void DelaunayGraphics::draw_vertices() {
		glDisable(GL_LIGHTING) ;
		glPointSize(GLfloat(vertices_size_ * 20)) ;
		int w,h ;
		glut_viewer_get_screen_size(&w, &h) ;
		notify_screen_resize(w,h) ;
		glEnable(GL_POINT_SMOOTH) ;
		glBegin(GL_POINTS) ;
		//        begin_spheres() ;
		if(!vertices_color_) {
			FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
				if(!show_copies_ && !delaunay_->is_primary(v))
					continue ;
				/*				if(v->locked) {
				glColor3f(1.0,0.0,0.0) ;
				} else*/ {
					glColor3f(0.0,0.0,0.0) ;
				}
				glVertex2f(v->point().x(), v->point().y()) ;
			}
		}
		else {
			FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
				if(!show_copies_ && !delaunay_->is_primary(v))
					continue ;
				int deg = delaunay_->degree(v) ;
				gl_vertex_color(deg) ;
				glVertex2f(v->point().x(), v->point().y()) ;
			}
		}
		//        end_spheres() ;
		glEnd() ;
		if(vertices_size_ > 0.6) {
			glBegin(GL_LINES) ;
			FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
				double R = 0.5 * radius(v) ;
				vec2 p = to_geex(v->point()) ;

				if(!show_copies_ && !delaunay_->is_primary(v))
					continue ;

				vec2 U,V ;
				CVT_->query_anisotropy(p,U,V) ;
				/*
				vec2 U(cos(v->theta), sin(v->theta)) ;
				U *= R ;
				vec2 V(-U.y, U.x) ; 
				*/              

				U *= 0.3 ;
				V *= 0.3 ;
				glVertex(p-U) ; glVertex(p+U) ;
				glVertex(p-V) ; glVertex(p+V) ;
			}
			glEnd() ;
		}
	}

	void DelaunayGraphics::draw_centers() {
		glDisable(GL_LIGHTING) ;
		glPointSize(int(centers_size_ * 20)) ;
		glColor3f(0.0,0.5,0.0) ;
		int w,h ;
		glut_viewer_get_screen_size(&w, &h) ;
		notify_screen_resize(w,h) ;
		glEnable(GL_POINT_SMOOTH) ;
		glBegin(GL_POINTS) ;
		//        begin_spheres() ;
		FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
			double V ;
			vec2 g ;

			if(!show_copies_ && !delaunay_->is_primary(v))
				continue ;

			if(delaunay_->period()) {
				if(!delaunay_->is_primary(v))
					continue ;
				CVT_->get_cell_primary_centroid(v, g, V) ;
			}
			else 
				CVT_->get_cell_centroid(v, g, V) ;
			vec2 p = g ;
			glVertex2f(p.x, p.y) ;
		}
		glEnd() ;
		//        end_spheres() ;
	}

	void DelaunayGraphics::draw_cells(bool colorize, bool mesh) {
		FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
			if(colorize) {gl_random_color() ; }
			if(delaunay_->dual_cell_intersects_boundary(v) || delaunay_->dimension()==1) { 
				Polygon2* P = delaunay_->dual_convex_clip(v) ;
				if(mesh) {
					glBegin(GL_LINES) ;
					for(unsigned int i=0; i<P->size(); i++) {
						glVertex((*P)[i].vertex[0]) ;
						glVertex((*P)[i].vertex[1]) ;
					}
					glEnd() ;
				} else {
					glBegin(GL_TRIANGLES) ;
					for(unsigned int i=0; i<P->size(); i++) {
						glVertex(to_geex(v->point())) ;
						glVertex((*P)[i].vertex[0]) ;
						glVertex((*P)[i].vertex[1]) ;
					}
					glEnd() ;
				}
			} else {
				glBegin(GL_POLYGON) ;
				Delaunay::Face_circulator it = delaunay_->incident_faces(v) ;
				do {
					glVertex(it->dual) ;
					it++ ;
				} while(it != delaunay_->incident_faces(v)) ;
				glEnd() ;
			} 
			//if(!delaunay_->dual_cell_infinite(v)) { 
			//	glBegin(GL_POLYGON) ;
			//	Delaunay::Face_circulator it = delaunay_->incident_faces(v) ;
			//	do {
			//		glVertex(it->dual) ;
			//		it++ ;
			//	} while(it != delaunay_->incident_faces(v)) ;
			//	glEnd() ;
			//}
		}
	}

	void DelaunayGraphics::draw_inner_voronoi() {
		glColor3f(0.0f, 1.0f, 0.0f) ;
		glLineWidth(2.0f) ;	
		glBegin(GL_LINES) ;
		FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
			if(!delaunay_->is_primary(v)) continue ;
			std::vector<vec2> P ;
			delaunay_->compute_inner_voronoi(v, P) ;
			if(P.size()<4) continue ;
			for(unsigned int i=0; i<P.size(); i++) {
				glVertex(P[i]) ;
				glVertex(P[(i+1)%P.size()]) ;
			}
		}
		glEnd() ;

		glPointSize(vertices_size_*20.f) ;
		glColor3f(1.0, 0.5, 0.0) ;
		glBegin(GL_POINTS) ;
		FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
			if(!delaunay_->is_primary(v)) continue ;
			std::vector<vec2> P ;
			delaunay_->compute_inner_voronoi(v, P) ;
			glVertex(polygon_centroid(P)) ;
		}
		glEnd() ;
	}

	void DelaunayGraphics::draw_dual_edges() {
		glColor3f(0.0f, 0.0f, 0.0f) ;
		glLineWidth(1.0f) ;

		FOR_EACH_EDGE_DT(Delaunay, delaunay_, e) {
			CGAL::Object o = delaunay_->baseclass::dual(e) ;
			Ray ray ;
			CGAL_Segment seg ;
			CGAL_Line line ;
			glBegin(GL_LINES) ;
			if(CGAL::assign(ray, o)) {
				glVertex(to_geex(ray.source())) ;
				glVertex(to_geex(ray.source()+1e3*ray.to_vector())) ;
			}
			else if(CGAL::assign(seg, o)) {
				glVertex(to_geex(seg.start())) ;
				glVertex(to_geex(seg.end())) ;
			}
			else if(CGAL::assign(line, o)) {
				if(line.b()!=0) {
					glVertex(vec2(-1e3, -1e3*(-line.a()/line.b())-line.c()/line.b())) ;
					glVertex(vec2(1e3,  1e3*(-line.a()/line.b())-line.c()/line.b())) ;
				} else {
					glVertex(vec2(-1e3*(-line.b()/line.a())-line.c()/line.a(), -1e3)) ;
					glVertex(vec2( 1e3*(-line.b()/line.a())-line.c()/line.a(),  1e3)) ;
				}
			}
			glEnd() ;
		}
	}

	void DelaunayGraphics::draw_non_hex_cells() {
		FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
			if(delaunay_->dual_cell_intersects_boundary(v)) { 
				Polygon2* P = delaunay_->dual_convex_clip(v) ;
				glColor3f(0.0, 0.5, 0.0) ;
				glBegin(GL_POLYGON) ;
				for(unsigned int i=0; i<P->size(); i++) {
					//                    glVertex(to_geex(v->point())) ;
					glVertex((*P)[i].vertex[0]) ;
					glVertex((*P)[i].vertex[1]) ;
				}
				glEnd() ;
			} else  {
				int degree = delaunay_->dual_facet_degree(v, delaunay_->period()) ;
				if(delaunay_->dual_facet_degree(v, delaunay_->period()) != 6 && !delaunay_->dual_cell_infinite(v)) {
					if(degree > 7)
						glColor3f(0.0, 0.0, 0.5) ;
					else if(degree == 7)
						glColor3f(0.0, 0.0, 1) ;
					else if(degree == 5)
						glColor3f(0.0, 0.64, 0.9) ;
					else
						glColor3f(0.0, 0.9, 0.9) ;
					glBegin(GL_POLYGON) ;
					Delaunay::Face_circulator it = delaunay_->incident_faces(v) ;
					do {
						glVertex(it->dual) ;
						it++ ;
					} while(it != delaunay_->incident_faces(v)) ;
					glEnd() ;
				}
			} 
		}
	}

	void DelaunayGraphics::draw_non_hex_periodic_cells() {
		FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
			int degree = delaunay_->dual_facet_degree_period(v) ;
			if( degree != 6) {
				if(degree > 7)
					glColor3f(0.0, 0.0, 0.5) ;
				else if(degree == 7)
					glColor3f(0.0, 0.0, 1) ;
				else if(degree == 5)
					glColor3f(0.0, 0.64, 0.9) ;
				else
					glColor3f(0.0, 0.9, 0.9) ;

				if(delaunay_->dual_cell_intersects_boundary(v) || delaunay_->dimension()==1) { 
					Polygon2* P = delaunay_->dual_convex_clip(v) ;
					glBegin(GL_POLYGON) ;
					for(unsigned int i=0; i<P->size(); i++) {
						glVertex((*P)[i].vertex[0]) ;
						glVertex((*P)[i].vertex[1]) ;
					}
					glEnd() ;
				}
				else if(delaunay_->is_primary(v)) {
					glBegin(GL_POLYGON) ;
					Delaunay::Face_circulator it = delaunay_->incident_faces(v) ;
					if(delaunay_->dimension()>1) {
						do {
							glVertex(it->dual) ;
							it++ ;
						} while(it != delaunay_->incident_faces(v)) ;
					}
					glEnd() ;
				}
			}
		}
	}

	inline bool is_quad(Delaunay* del, Delaunay::Face_handle f, int i, double ratio) {
		Delaunay::Face_handle g = f->neighbor(i) ;
		if(del->is_infinite(g)) { return false ; }

		vec2 p1 = to_geex(f->vertex(0)->point()) ;
		vec2 q1 = to_geex(g->vertex(0)->point()) ;

		vec2 c1 = f->dual ;
		vec2 c2 = g->dual ;

		double R1 = (p1 - c1).length() ;
		double R2 = (q1 - c2).length() ;
		double R = 0.5 * (R1 + R2) ;

		double dC = (c2-c1).length() ;
		return (dC < ratio * R) ;
	}

	bool is_obtuse(Delaunay::Face_handle fh) {
		vec2 v1 = to_geex(fh->vertex(0)->point()) ;
		vec2 v2 = to_geex(fh->vertex(1)->point()) ;
		vec2 v3 = to_geex(fh->vertex(2)->point()) ;
		double d1 = distance2(v2, v3) ;
		double d2 = distance2(v3, v1) ;
		double d3 = distance2(v1, v2) ;
		return d1+d2<d3 || d2+d3<d1 || d3+d1<d2 ;
	}

	void DelaunayGraphics::draw_primal_mesh() {
		glLineWidth(1.0) ;
		glColor3f(0.0, 0.0, 0.9) ;
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE) ;
		glBegin(GL_LINES) ;
		FOR_EACH_FACE_DT(Delaunay, delaunay_, it) {
			if(delaunay_->is_infinite(it)) { continue ; }
			//	if(it->dual_outside) continue ;

			/* we don't draw quad here */
			for(unsigned int i=0; i<3; i++) {
				unsigned int j1 = i + 1 ; 
				if(j1 == 3) { j1 = 0 ; }
				unsigned int j2 = (j1 + 1) ;
				if(j2 == 3) { j2 = 0 ; }

				//if(!is_quad(delaunay_, it, i, quad_ratio_)) {
				Delaunay::Vertex_handle v1 = it->vertex(j1) ;
				Delaunay::Vertex_handle v2 = it->vertex(j2) ;
				glVertex2f(v1->point().x(), v1->point().y()) ;
				glVertex2f(v2->point().x(), v2->point().y()) ;
				//}
			}

		}
		glEnd() ;
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL) ;
		glColor3f(0.8, 0.8, 0.8) ;
		glBegin(GL_TRIANGLES) ;
		FOR_EACH_FACE_DT(Delaunay, delaunay_, it) {
			if(delaunay_->is_infinite(it)) { continue ; }
			if(is_obtuse(it)) {
				for(int i=0; i<3; ++i) {
					Delaunay::Vertex_handle vi = it->vertex(i) ;
					glVertex(to_geex(vi->point())) ;
				}
			}
		}
		glEnd() ;
	}

	void DelaunayGraphics::draw_dual_mesh() {

		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE) ;
		//glLineWidth(3) ;
		//glColor3f(0.0, 0.0, 0.0) ;
		//draw_cells(false, true) ;
		//glLineWidth(9) ;
		//glColor3f(1.0, 1.0, 1.0) ;
		glColor3f(0.0, 0.0, 0.0) ;
		glLineWidth(1) ;
		draw_cells(false, true) ;
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL) ;
	}

	void DelaunayGraphics::is_min_max(Delaunay::Vertex_handle v, bool& is_min, bool& is_max) {
		is_min = true ;
		is_max = true ;
		Delaunay::Vertex_circulator it = delaunay_->incident_vertices(v) ;
		do {
			if(!delaunay_->is_infinite(it)) {
				//                is_min = is_min && v->energy < it->energy ;
				//                is_max = is_max && v->energy > it->energy ;
			}
			it++ ;
		} while(it != delaunay_->incident_vertices(v)) ;
	}

	///dxy change: draw energy (hsv, green to red)
	void DelaunayGraphics::draw_energy() {

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL) ;

		double emin = 1e30 ;
		double emax = -1e30 ;

		Delaunay::Vertex_handle vmin = 0 ;
		Delaunay::Vertex_handle vmax = 0 ;

		// 		FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) { //dxy change
		FOR_EACH_PRIMARY_VERTEX_DT(Delaunay, delaunay_, v) {
			if(v->energy < emin) { emin = v->energy; vmin = v ; }
			if(v->energy > emax) { emax = v->energy; vmax = v ; }
		}

		//        std::cerr << "emin = " << emin << " emax =" << emax << std::endl ;

		double emean = 0.5 * (emin + emax) ;


		/* if(emax - emean > emean - emin) */ {
			//glColor3f(1.0, 0.5, 0.5) ;
			glColor3f(1.0, 0.0, 0.0) ;
			if(vmax != 0) { /*draw_dual_facet(vmax)*/ draw_period_dual_facet(vmax) ; }
		} 
		/* else */ {
			//glColor3f(0.5, 0.5, 1.0) ;
			glColor3f(0.0, 1.0, 0.0);
			if(vmin != 0) { /*draw_dual_facet(vmin)*/ draw_period_dual_facet(vmin) ; }
		}


		double scale = (emax - emin) ;
		if(::fabs(scale) < 1e-30) { scale = 1.0 ; }
		scale = 1.0 / scale ;

		///dxy change
		const vec3 c1(120, 1.0, 1.0);
		const vec3 c2(0.0, 1.0, 1.0);
		///
		//         const vec3 c1(0.0, 0.0, 1.0) ;
		//         const vec3 c2(1.0, 0.0, 0.0) ;

		FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
			double e = scale * (v->energy - emin) ;
			//glColor(mix(c1, c2, e)) ;
			///dxy add
			vec3 mix_hsv = mix(c1, c2, e);
			vec3 mix_rgb;
			hsv2rgb(mix_hsv, mix_rgb);
			glColor3f(mix_rgb.x, mix_rgb.y, mix_rgb.z);
			///
			/*draw_dual_facet(v) ;*/
			draw_period_dual_facet(v);
		}
	}

	void DelaunayGraphics::draw_regularity() {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL) ;

		double emin = 1e30 ;
		double emax = -1e30 ;

		Delaunay::Vertex_handle vmin = 0 ;
		Delaunay::Vertex_handle vmax = 0 ;

		FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
			if(!delaunay_->is_primary(v))
				continue ;
			if(v->regularity < emin) { emin = v->regularity; vmin = v ; }
			if(v->regularity > emax) { emax = v->regularity; vmax = v ; }
		}

		//        std::cerr << "emin = " << emin << " emax =" << emax << std::endl ;

		double emean = 0.5 * (emin + emax) ;

		// 2d optimital measure
		emin = 5.0/(36*sqrt(3.0)) ;
		//		emax = 5.0/(36*sqrt(3.0)) * 1.2;

		///* if(emax - emean > emean - emin) */ {
		//    glColor3f(1.0, 0.5, 0.5) ;
		//    if(vmax != 0) { draw_dual_facet(vmax) ; }
		//} 
		///* else */ {
		//    glColor3f(0.5, 0.5, 1.0) ;
		//    if(vmin != 0) { draw_dual_facet(vmin) ; }
		//}

		double scale = (emax - emin) ;
		if(::fabs(scale) < 1e-30) { scale = 1.0 ; }
		scale = 1.0 / scale ;

		//		std::cout << "max/min ratio: " << emax/emin << std::endl ;

		const vec3 c1(0.0, 0.0, 1.0) ;
		const vec3 c2(1.0, 0.0, 0.0) ;

		FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
			double e ;
			if(delaunay_->is_primary(v)) {
				e = scale * (v->regularity - emin) ;
			}
			else {
				Delaunay::Vertex_handle vp = delaunay_->vertices()[v->index] ;
				e = scale * (vp->regularity - emin) ;
			}
			glColor(mix(c1, c2, e)) ;
			//            draw_dual_facet(v) ;

			if(delaunay_->dual_cell_intersects_boundary(v) || delaunay_->dimension()==1) { 
				Polygon2* P = delaunay_->dual_convex_clip(v) ;
				glBegin(GL_POLYGON) ;
				for(unsigned int i=0; i<P->size(); i++) {
					glVertex((*P)[i].vertex[0]) ;
					glVertex((*P)[i].vertex[1]) ;
				}
				glEnd() ;
			}
			else if(delaunay_->is_primary(v)) {
				glBegin(GL_POLYGON) ;
				Delaunay::Face_circulator it = delaunay_->incident_faces(v) ;
				if(delaunay_->dimension()>1) {
					do {
						glVertex(it->dual) ;
						it++ ;
					} while(it != delaunay_->incident_faces(v)) ;
				}
				glEnd() ;
			}
		}
	}

	void DelaunayGraphics::draw_field() {
		double x_min, y_min, z_min ;
		double x_max, y_max, z_max ;
		delaunay_->get_bbox(x_min, y_min, z_min, x_max, y_max, z_max) ;
		double dx = x_max - x_min ;
		double dy = y_max - y_min ;
		glPointSize(4) ;
		glBegin(GL_POINTS) ;

		for(int i=0; i<100; i++) {
			for(int j=0; j<100; j++) {
				vec2 p(x_min + dx * double(i)/100.0, y_min + dy * double(j)/100.0) ;
				gl_random_color(delaunay_->segments().locate(p)) ;
				glVertex(p) ;
			}
		}

		glPointSize(8) ;
		for(SegmentDelaunay::Vertex_iterator it = delaunay_->segments().finite_vertices_begin() ;
			it != delaunay_->segments().finite_vertices_end(); it++
			) {
				gl_random_color(it->index) ;
				glVertex(to_geex(it->point())) ;
		}

		glEnd() ;
	}

	void DelaunayGraphics::draw_dual_facet(Delaunay::Vertex_handle v) {
		if(delaunay_->dual_cell_intersects_boundary(v)) { 
			Polygon2* P = delaunay_->dual_convex_clip(v) ;
			glBegin(GL_TRIANGLES) ;
			for(unsigned int i=0; i<P->size(); i++) {
				glVertex(to_geex(v->point())) ;
				glVertex((*P)[i].vertex[0]) ;
				glVertex((*P)[i].vertex[1]) ;
			}
			glEnd() ;
		} else  {
			glBegin(GL_POLYGON) ;
			Delaunay::Face_circulator it = delaunay_->incident_faces(v) ;
			do {
				glVertex(it->dual) ;
				it++ ;
			} while(it != delaunay_->incident_faces(v)) ;
			glEnd() ;
		} 
	}

	///dxy add
	void DelaunayGraphics::draw_period_dual_facet(Delaunay::Vertex_handle v) {
		if (!delaunay_->period())
		{
			draw_dual_facet(v);
			return;
		}

		if (!delaunay_->is_primary(v)) return;
		if(delaunay_->dual_cell_intersects_boundary(v)) { 
			Polygon2* P = delaunay_->dual_convex_clip(v) ;
			glBegin(GL_POLYGON);
			for(unsigned int i=0; i<P->size(); i++) {
				//const vec2& p1 = (*P)[i].vertex[0] ;
				//const vec2& p2 = (*P)[i].vertex[1] ;
				glVertex((*P)[i].vertex[0]);
				glVertex((*P)[i].vertex[1]);
			}
			glEnd();

			//mirror sites.
			std::set<Delaunay::Vertex_handle>& m = delaunay_->mirrors(v->index) ;
			for(std::set<Delaunay::Vertex_handle>::iterator it=m.begin(); it!=m.end(); ++it) {
				Polygon2* P = delaunay_->dual_convex_clip(*it) ;
				glBegin(GL_POLYGON);
				for (unsigned int i=0; i<P->size(); ++i)
				{
					glVertex((*P)[i].vertex[0]);
					glVertex((*P)[i].vertex[1]);
				}
				glEnd();
				/*vec2 offset = p0 - to_geex((*it)->point()) ;
				for(unsigned int i=0; i<P->size(); ++i) {
				const vec2& p1 = (*P)[i].vertex[0] + offset ;
				const vec2& p2 = (*P)[i].vertex[1] + offset ;
				double Vi = triangle_area(p0, p1, p2) ;
				vec2 Gi = triangle_centroid(p0, p1, p2) ;
				V += Vi ;
				Vg += Vi * Gi ;
				lloyd += Lloyd_energy(p0, p1, p2) ;
				}*/
			}
		} else {
			Delaunay::Face_circulator it = delaunay_->incident_faces(v) ;
			glBegin(GL_POLYGON);
			do 
			{
				glVertex(it->dual);
				it++;
			} while (it != delaunay_->incident_faces(v));
			glEnd();
			/*Delaunay::Face_circulator jt = it ; jt++ ;
			do {
			const vec2& p1 = it->dual ;
			const vec2& p2 = jt->dual ;
			double Vi = triangle_area(p0, p1, p2) ;
			vec2 Gi = triangle_centroid(p0, p1, p2) ;
			V += Vi ;
			Vg += Vi * Gi ;
			lloyd += Lloyd_energy(p0, p1, p2) ;
			it++ ; jt++ ;
			} while(it != delaunay_->incident_faces(v)) ;*/
		}
	}
	///


	double DelaunayGraphics::radius(Delaunay::Vertex_handle v) {
		double result = 1e30 ;
		vec2 p1 = to_geex(v->point()) ;
		Delaunay::Vertex_circulator it = delaunay_->incident_vertices(v) ;
		do {
			if(!delaunay_->is_infinite(it)) {
				vec2 p2 = to_geex(it->point()) ;
				result = gx_min(result, (p2 - p1).length2()) ;
			}
			it++ ;
		} while(it != delaunay_->incident_vertices(v)) ;
		return sqrt(result) ;
	}

	void DelaunayGraphics::draw_boundary_cells() {
		int current_color = 0 ;
		//		delaunay_->compute_rvd() ;

		FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
			if(!delaunay_->is_primary(v)) continue ;

			if(v->dual_intersects_boundary || delaunay_->dimension()==1) {
				if(colorize_) {
					current_color = random_color_index_ ;
					gl_random_color() ;
				}
				else glColor3f(1.f, 1.f, 0.0f) ;

				Polygon2 *pl = delaunay_->dual_convex_clip(v) ; //rvd_[v->index] ;
				glBegin(GL_POLYGON) ; // convex only, need to improve
				for(unsigned int j=0; j<pl->size(); ++j) {
					glVertex((*pl)[j].vertex[0]) ;
					glVertex((*pl)[j].vertex[1]) ;
				}
				glEnd() ;

				// draw border
				glColor3f(0.f, 0.f, 0.f) ;
				glBegin(GL_LINES) ;
				for(unsigned int j=0; j<pl->size(); ++j) {
					glVertex((*pl)[j].vertex[0]) ;
					glVertex((*pl)[j].vertex[1]) ;
				}
				glEnd() ;

				// draw mirrors
				std::set<Delaunay::Vertex_handle>& mirrors = delaunay_->mirrors(v->index) ;
				for(std::set<Delaunay::Vertex_handle>::iterator it=mirrors.begin(); it!=mirrors.end(); ++it) {
					if((*it)->dual_intersects_boundary) {
						Polygon2 *pl = delaunay_->dual_convex_clip(*it) ; //rvd_[v->index] ;
						if(colorize_) gl_random_color(current_color) ;

						if(show_pvd_euclidean_) { // offset polygon
							//int mid = delaunay_->domain_idx(to_geex((*it)->point())) ;
							vec2 offset = to_geex(v->point()) - to_geex((*it)->point());
							for(unsigned int i=0; i<pl->size(); ++i) {
								(*pl)[i].vertex[0] = (*pl)[i].vertex[0] + offset ;
								(*pl)[i].vertex[1] = (*pl)[i].vertex[1] + offset ;
							}
						}

						glBegin(GL_POLYGON) ; // convex only, need to improve
						for(unsigned int j=0; j<pl->size(); ++j) {
							glVertex((*pl)[j].vertex[0]) ;
							glVertex((*pl)[j].vertex[1]) ;
						}
						glEnd() ;

						// draw border
						glColor3f(0.f, 0.f, 0.f) ;
						glBegin(GL_LINES) ;
						for(unsigned int j=0; j<pl->size(); ++j) {
							glVertex((*pl)[j].vertex[0]) ;
							glVertex((*pl)[j].vertex[1]) ;
						}
						glEnd() ;
					}
				}
			}
			else {
				glColor3f(0.0, 0.0, 0.0) ;
				glBegin(GL_LINE_LOOP) ;
				Delaunay::Face_circulator it = delaunay_->incident_faces(v) ;
				if(delaunay_->dimension()>1) { // to avoid degenerate case: all the points are on a line.
					do {
						glVertex(it->dual) ;
						it++ ;
					} while(it != delaunay_->incident_faces(v)) ;
				}
				glEnd() ;
			}
		}
	}

	void DelaunayGraphics::draw_edge_hist() {
		std::vector<double>& hist = IO_->edge_histogram() ;
		int w, h ;
		glut_viewer_get_screen_size(&w, &h) ;
		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();
		gluOrtho2D(0, w, 0, h);

		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glLoadIdentity();

		double maxhist = 0 ;
		for(unsigned int i=0; i<hist.size(); ++i) {
			if(maxhist < hist[i])
				maxhist = hist[i] ;
		}

		glColor4f(0.76, 0.87, 0.85, 0.5) ;
		glBegin(GL_QUADS) ;
		glVertex2i(w-256, h-256) ;
		glVertex2i(w, h-256) ;
		glVertex2i(w, h) ;
		glVertex2i(w-256, h) ;
		glEnd() ;

		glLineWidth(3.0) ;
		glColor3f(0.2, 0.2, 0.2) ;
		glBegin(GL_LINES) ;
		for(int i=0; i<256; ++i) {
			glVertex2i(w-256+i,   h-256) ;
			glVertex2i(w-256+i,   h-256+hist[i]/maxhist*256) ;
		}
		glEnd() ;

		glMatrixMode(GL_PROJECTION);
		glPopMatrix() ;

		glMatrixMode(GL_MODELVIEW);
		glPopMatrix() ;
	}

	static void drawCircle(float cx, float cy, float r, int num_segments) 
	{ 
		glBegin(GL_LINE_LOOP); 
		for(int ii = 0; ii < num_segments; ii++) 
		{ 
			float theta = 2.0f * 3.1415926f * float(ii) / float(num_segments);//get the current angle 
			float x = r * cosf(theta);//calculate the x component 
			float y = r * sinf(theta);//calculate the y component 
			glVertex2f(x + cx, y + cy);//output vertex 
		} 
		glEnd(); 
	}

	static inline void draw_vertex_disk(GLUquadric *q, const Delaunay::Vertex_handle& v, double radius) {
		glPushMatrix() ;
		glTranslatef(v->point().x(), v->point().y(), 0) ;
		glColor3f(.7, 0.7, 0.7) ;
		gluDisk(q, 0, radius, 100, 100) ;
		glPopMatrix() ;
		glColor3f(0.0, 0.0, 0.0) ;
		drawCircle(v->point().x(), v->point().y(), radius, 100) ;
	}

	void DelaunayGraphics::draw_disk() {
		GLUquadric * q = gluNewQuadric( ) ;
		glLineWidth(2.0) ;

		FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
			if(!delaunay_->is_primary(v))
				continue ;
			draw_vertex_disk(q, v, delaunay_->sample_radius()) ;
		}
		gluDeleteQuadric(q) ;

		GLint func ;
		glGetIntegerv(GL_DEPTH_FUNC, &func) ;
		glDepthFunc(GL_LEQUAL) ;

		// draw skeleton of the gaps
		std::vector<std::vector<vec2> >& skeleton = delaunay_->void_skeleton() ;
		glLineWidth(4.0) ;
		glColor3f(0.0, 1.0, 0.0) ;
		for(int i=0; i<skeleton.size(); ++i) {
			glBegin(GL_LINES) ;
			for(unsigned int j=0; j<skeleton[i].size()-1; j+=2)  {
				glVertex(skeleton[i][j]) ;
				glVertex(skeleton[i][j+1]) ;
			}
			glEnd() ;
		}

		// draw gaps
		std::vector<std::vector<vec2> >& regions = delaunay_->void_regions() ;
		for(int i=0; i<regions.size(); ++i) {
			gl_random_color() ;
			glBegin(GL_POLYGON) ;
			for(unsigned int j=0; j<regions[i].size(); ++j)  {
				glVertex(regions[i][j]) ;
			}
			glEnd() ;
		}
		glDepthFunc(func) ;
	}

	void DelaunayGraphics::draw_min_max() {
		GLUquadric * q = gluNewQuadric( ) ;
		glLineWidth(2.0) ;

		// draw triangle with angle < 30
		FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
			if(!show_copies_ && !(delaunay_->is_primary(v) || delaunay_->neighbor_to_primary(v)))
				continue ;

			vec2 p = to_geex(v->point()) ;
			Delaunay::Vertex_circulator cir = delaunay_->incident_vertices(v) ;
			Delaunay::Vertex_circulator cir2 = cir++ ;
			do {
				vec2 e1 = to_geex(cir->point())-p ;
				vec2 e2 = to_geex(cir2->point())-p ;
				double angle = to_degree(acos(dot(e1, e2)/(e1.length()*e2.length()))) ;
				if(angle < 30) {
					glColor3f(1.0, 0, 0) ;
					glBegin(GL_POLYGON) ;
					glVertex(p) ;
					glVertex(to_geex(cir->point())) ;
					glVertex(to_geex(cir2->point())) ;
					glEnd() ;

					draw_vertex_disk(q, v, delaunay_->sample_radius()) ;
					draw_vertex_disk(q, cir, delaunay_->sample_radius()) ;
					draw_vertex_disk(q, cir2, delaunay_->sample_radius()) ;
				}
				++cir  ;
				++cir2 ;
			} while(cir!=delaunay_->incident_vertices(v)) ;
		}
		gluDeleteQuadric(q) ;

		// draw faces with max/min area
		Delaunay::Face_handle fmin, fmax ;
		double area_min=1e10, area_max=-1e10 ;

		FOR_EACH_FINITE_FACE_DT(Delaunay, delaunay_, f) {
			if(!f->dual_outside) {
				double area = triangle_area(to_geex(f->vertex(0)->point()), 
					to_geex(f->vertex(1)->point()), 
					to_geex(f->vertex(2)->point())) ;
				if(area < area_min) {
					area_min = area ;
					fmin = f ;
				}
				if(area > area_max) {
					area_max = area ;
					fmax = f ;
				}
			}
		}

		glColor3f(1.0, 0, 0) ;
		glBegin(GL_POLYGON) ;
		glVertex(to_geex(fmin->vertex(0)->point())) ;
		glVertex(to_geex(fmin->vertex(1)->point())) ;
		glVertex(to_geex(fmin->vertex(2)->point())) ;
		glEnd() ;

		glBegin(GL_POLYGON) ;
		glVertex(to_geex(fmax->vertex(0)->point())) ;
		glVertex(to_geex(fmax->vertex(1)->point())) ;
		glVertex(to_geex(fmax->vertex(2)->point())) ;
		glEnd() ;
	} 

	void DelaunayGraphics::draw_mesh() {
		Map* M = delaunay_->map() ;

		///dxy test
// 		std::cout << "Map: " << M->nb_vertices() << " vertices, " << M->nb_facets() << " facets." << std::endl;
// 		if (M->nb_vertices() == 0)
// 		{
// 			delaunay_->editor()->delaunay_to_map();
// 			std::cout << "Map: " << M->nb_vertices() << " vertices, " << M->nb_facets() << " facets." << std::endl;
// 		}
		///
		GLint func ;
		glGetIntegerv(GL_DEPTH_FUNC, &func) ;
		glDepthFunc(GL_LEQUAL) ;
		glLineWidth(2.0) ;
		glColor3f(0.2, 0.2, 0.2) ;
		glBegin(GL_LINES) ;
		for(Map::Facet_iterator it=M->facets_begin(); it!=M->facets_end(); ++it) {
			Map::Halfedge *h = it->halfedge() ;
			do {
				glVertex(h->prev()->vertex()->point()) ;
				glVertex(h->vertex()->point()) ;
				h = h->next() ;
			} while(h!=it->halfedge()) ;

		}
		glEnd() ;

		if(vertices_size_>0) {
			glPointSize(GLfloat(vertices_size_ * 30)) ;
			glBegin(GL_POINTS) ;
			for(Map::Vertex_iterator it=M->vertices_begin(); it!=M->vertices_end(); ++it) {
				gl_vertex_color(it->degree()) ;
				glVertex(it->point()) ;
			}
			glEnd() ;
		}

		// draw picked vertices and edges 
		IrregEditor *editor = delaunay_->editor() ;
		std::vector<Map::Vertex*>&   vsel = editor->picked_verts() ;
		std::vector<Map::Halfedge*>& esel = editor->picked_edges() ; 
		GLUquadric *q = gluNewQuadric() ;

		glEnable(GL_LIGHTING) ;
		glColor3f(0.1, 0.8, 0.8) ;
		for(unsigned int i=0; i<esel.size(); ++i) {
			Map::Halfedge * h = esel[i] ;
			vec3 v0 = h->prev()->vertex()->point() ;
			vec3 v1 = h->vertex()->point() ;
			vec3 dir = normalize(v1-v0) ;
			
			glPushMatrix() ;
			glTranslatef(v0.x, v0.y, 0) ;
			glRotatef(90, -dir.y, dir.x, 0.0) ;
			gluCylinder(q, 0.01, 0.01, distance(v0, v1), 20, 20) ;
			glPopMatrix() ;
		}

		glDisable(GL_LIGHTING) ;
		for(unsigned int i=0; i<vsel.size(); ++i) {
			vec3 c = vsel[i]->point() ;
			int  deg = vsel[i]->degree() ;
			
			gl_vertex_color(deg) ;
			glPushMatrix() ;
			glTranslatef(c.x, c.y, 0) ;
			gluSphere(q, vertices_size_/10.0, 20, 20) ;
			glPopMatrix() ;
		}

		gluDeleteQuadric(q) ;

		// draw current facet where the picked point locates
		Map::Facet* f = editor->cur_facet() ;
		if(f!=nil) {
			glColor3f(1.0, 0, 0) ;
			glBegin(GL_POLYGON) ;
			Map::Halfedge* h = f->halfedge() ;
			do{
				glVertex(h->vertex()->point()) ;
				h = h->next() ;
			} while(h!= editor->cur_facet()->halfedge()) ;
			glEnd() ;
		}

		glDepthFunc(func) ;
	}

	void DelaunayGraphics::draw_selection() {
		std::vector<Delaunay::Vertex_handle>&  vsel = delaunay_->vertices_sel() ;
		std::vector<Delaunay::Edge>& esel = delaunay_->edges_sel() ; 
		GLUquadric *q = gluNewQuadric() ;

		glEnable(GL_LIGHTING) ;
		glColor3f(0.1, 0.8, 0.8) ;
		for(unsigned int i=0; i<esel.size(); ++i) {
			Delaunay::Face_handle f = esel[i].first ;
			int  idx = esel[i].second ;
			vec2 v0 = to_geex(f->vertex(f->cw(idx))->point()) ;
			vec2 v1 = to_geex(f->vertex(f->ccw(idx))->point()) ;
			vec2 dir = normalize(v1-v0) ;

			glPushMatrix() ;
			glTranslatef(v0.x, v0.y, 0) ;
			glRotatef(90, -dir.y, dir.x, 0.0) ;
			gluCylinder(q, vertices_size_, vertices_size_, distance(v0, v1), 20, 20) ;
			glPopMatrix() ;
		}

		//		glDisable(GL_LIGHTING) ;
		for(unsigned int i=0; i<vsel.size(); ++i) {
			vec2 c = to_geex(vsel[i]->point()) ;
			int  deg = delaunay_->degree(vsel[i]) ;

			gl_vertex_color(deg) ;
			glPushMatrix() ;
			glTranslatef(c.x, c.y, 0) ;
			gluSphere(q, vertices_size_, 20, 20) ;
			glPopMatrix() ;
		}

		gluDeleteQuadric(q) ;
		//dxy add
		glDisable(GL_LIGHTING);
	}

	//////////////////////////////////////////////////////////////////////////
	///dxy add
	void DelaunayGraphics::draw_constrained_edge() {
		glLineWidth(1) ;
		glColor3f(0.0, 0.0, 0.5) ;
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE) ;
		glBegin(GL_LINES) ;
		FOR_EACH_FACE_DT(Delaunay, delaunay_, it) {
			if(delaunay_->is_infinite(it)) { continue ; }

			for(unsigned int i=0; i<3; i++) {
				unsigned int j1 = i + 1 ; 
				if(j1 == 3) { j1 = 0 ; }
				unsigned int j2 = (j1 + 1) ;
				if(j2 == 3) { j2 = 0 ; }

				if(!is_quad(delaunay_, it, i, quad_ratio_)) {
					Delaunay::Vertex_handle v1 = it->vertex(j1) ;
					Delaunay::Vertex_handle v2 = it->vertex(j2) ;
					int cw_j1 = (j1 + 2) % 3;
					Delaunay::Face_handle f = it->neighbor(cw_j1);
					if (delaunay_->is_infinite(f))
					{
						glColor3f(0.0, 0.0, 0.0);   //boundary edge: black
						glVertex2f(v1->point().x(), v1->point().y()) ;
						glVertex2f(v2->point().x(), v2->point().y()) ;
					}
					else {
						vec2 p1 = it->dual;
						vec2 p2 = f->dual;
						if (!delaunay_->in_boundary(p1) || !delaunay_->in_boundary(p2))
						{
							glColor3f(1.0, 0.0, 0); //constrained edge: red
							glVertex2f(v1->point().x(), v1->point().y()) ;
							glVertex2f(v2->point().x(), v2->point().y()) ;
						}
					}
				}
			}
		}
		glEnd() ;
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL) ;
	}

	void DelaunayGraphics::draw_lloyd_grad() {
		glLineWidth(3) ;
		glColor3f(0.0, 0.0, 1.0) ;
		FOR_EACH_VERTEX_DT(Delaunay, delaunay_, it) {
			if(delaunay_->is_infinite(it)) {continue ; }
			vec2 gstart = to_geex(it->point()) ;
			vec2 grad = -10000 * lloyd_grad_magnification_ * it->lloyd_grad;
			draw_vector(grad, gstart);
		}
	}

	void DelaunayGraphics::draw_direction_grad() {
		glLineWidth(3) ;
		glColor3f(1.0, 0.0, 0.0) ;
		glBegin(GL_LINES) ;
		FOR_EACH_VERTEX_DT(Delaunay, delaunay_, it) {
			if(delaunay_->is_infinite(it)) {continue ; }
			vec2 gstart = to_geex(it->point()) ;
			vec2 grad = -10000 * direct_grad_magnification_ * it->direction_grad ;
			draw_vector(grad, gstart);
		}
		glEnd() ;
	}

	void DelaunayGraphics::draw_regular_direction_grad() {
		delaunay_->calc_vertex_group_mark();

		glLineWidth(3) ;
		glColor3f(1.0, 0.0, 0.0) ;
		glBegin(GL_LINES) ;
		FOR_EACH_VERTEX_DT(Delaunay, delaunay_, it) {
			if (it->group_mark != -1) continue; //it belongs to some irregular face group
			
			if(delaunay_->is_infinite(it)) continue ;
			vec2 gstart = to_geex(it->point()) ;
			vec2 grad = -10000 * direct_grad_magnification_ * it->direction_grad ;
			draw_vector(grad, gstart);
		}
		glEnd() ;
	}

	void DelaunayGraphics::draw_total_grad() {
		glLineWidth(3) ;
		glColor3f(0.0, 1.0, 0.0) ;
		glBegin(GL_LINES) ;
		FOR_EACH_VERTEX_DT(Delaunay, delaunay_, it) {
			if(delaunay_->is_infinite(it)) {continue ; }
			vec2 gstart = to_geex(it->point()) ;
			vec2 grad = -10000 * total_grad_magnification_ * (it->direction_grad + it->lloyd_grad) ;
			draw_vector(grad, gstart);
		}
		glEnd() ;
	}

	void DelaunayGraphics::draw_direction_field() {
		double average_length = sqrt(delaunay_->average_area_);
		glLineWidth(1.5);
		glColor3f(0.0, 0.0, 0.0);
		FOR_EACH_VERTEX_DT(Delaunay, delaunay_, it) {
			vec2 start = to_geex(it->point());
			double field_angle;
			CVT_->get_field(start, field_angle);
			vec2 field(cos(field_angle), sin(field_angle));
			field *= 0.7 * average_length;
			draw_vector(field, start);
		}
	}

	//
	void DelaunayGraphics::draw_relative_area() {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL) ;

		double emin = 1e30 ;
		double emax = -1e30 ;

		Delaunay::Vertex_handle vmin = 0 ;
		Delaunay::Vertex_handle vmax = 0 ;

		// 		FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v) {
		FOR_EACH_PRIMARY_VERTEX_DT(Delaunay, delaunay_, v) {
			{
				if(v->area < emin) { emin = v->area; vmin = v ; }
				if(v->area > emax) { emax = v->area; vmax = v ; }
			}
		}

		double emean = 0.5 * (emin + emax) ;

		/* if(emax - emean > emean - emin) */ {
			glColor3f(1.0, 0.0, 1.0);
			if(vmax != 0) { draw_period_dual_facet(vmax)/*draw_dual_facet(vmax)*/ ; }
		} 
		/* else */ {
			glColor3f(0.0, 1.0, 1.0);
			if(vmin != 0) {/* draw_dual_facet(vmin)*/ draw_period_dual_facet(vmin); }
		}

		double scale = (emax - emin) ;
		if(::fabs(scale) < 1e-30) { scale = 1.0 ; }
		scale = 1.0 / scale ;

		const vec3 c1(120, 1.0, 1.0);
		const vec3 c2(0.0, 1.0, 1.0);

		/*FOR_EACH_VERTEX_DT(Delaunay, delaunay_, v)*/
		FOR_EACH_PRIMARY_VERTEX_DT(Delaunay, delaunay_, v)
		{
			double e = scale * (v->area - emin) ;
			vec3 mix_hsv = mix(c1, c2, e);
			vec3 mix_rgb;
			hsv2rgb(mix_hsv, mix_rgb);
			glColor3f(mix_rgb.x, mix_rgb.y, mix_rgb.z);
			/*draw_dual_facet(v) ;*/
			draw_period_dual_facet(v);
		}
	}

	void DelaunayGraphics::draw_primal_area() {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL) ;

		double emin = 1e30 ;
		double emax = -1e30 ;

		Delaunay::Face_handle fmin = 0 ;
		Delaunay::Face_handle fmax = 0 ;

		FOR_EACH_FACE_DT(Delaunay, delaunay_, f) {
			if (delaunay_->is_infinite(f)) continue;

			double farea = area(f->vertex(0)->point(),
				f->vertex(1)->point(),
				f->vertex(2)->point()
				);
			if(farea < emin) { emin = farea; fmin = f ; }
			if(farea > emax) { emax = farea; fmax = f ; }
		}
		double emean = 0.5 * (emin + emax) ;

		/* if(emax - emean > emean - emin) */ {
			glColor3f(1.0, 0.0, 1.0);
			if(fmax != 0) { draw_primal_triangle(fmax) ; }
		} 
		/* else */ {
			glColor3f(0.0, 1.0, 1.0);
			if(fmin != 0) { draw_primal_triangle(fmin) ; }
		}

		double scale = (emax - emin) ;
		if(::fabs(scale) < 1e-30) { scale = 1.0 ; }
		scale = 1.0 / scale ;

		const vec3 c1(120, 1.0, 1.0);
		const vec3 c2(0.0, 1.0, 1.0);

		FOR_EACH_FACE_DT(Delaunay, delaunay_, f) {
			double farea = area(f->vertex(0)->point(),
				f->vertex(1)->point(),
				f->vertex(2)->point()
				);
			double e = scale * (farea - emin) ;
			vec3 mix_hsv = mix(c1, c2, e);
			vec3 mix_rgb;
			hsv2rgb(mix_hsv, mix_rgb);
			glColor3f(mix_rgb.x, mix_rgb.y, mix_rgb.z);
			draw_primal_triangle(f) ;
		}
	}

	//
	void DelaunayGraphics::draw_long_edge()
	{
		int count = 0;
		glLineWidth(3) ;
		glColor3f(0.0, 0.0, 0.5) ;
		glBegin(GL_LINES) ;		
		FOR_EACH_FACE_DT(Delaunay, delaunay_, f) {
			if (delaunay_->is_infinite(f)) continue;
			for (int i=0; i<3; ++i)
			{
				Delaunay::Vertex_handle v1 = f->vertex((i+1)%3);
				Delaunay::Vertex_handle v2 = f->vertex((i+2)%3);
				if (delaunay_->is_long_edge(v1, v2))
				{
					count++;
					glVertex(to_geex(v1->point()));
					glVertex(to_geex(v2->point()));
				}
			}
		}
		glEnd();
		std::cout << count/2 << " long edges." << std::endl;
	}

	void DelaunayGraphics::draw_short_edge()
	{
		int count = 0;
		glLineWidth(3) ;
		glColor3f(0.5, 0.0, 0.0) ;
		glBegin(GL_LINES) ;		
		FOR_EACH_FACE_DT(Delaunay, delaunay_, f) {
			if (delaunay_->is_infinite(f)) continue;
			for (int i=0; i<3; ++i)
			{
				Delaunay::Vertex_handle v1 = f->vertex((i+1)%3);
				Delaunay::Vertex_handle v2 = f->vertex((i+2)%3);
				if (delaunay_->is_short_edge(v1, v2))
				{
					count++;
					glVertex(to_geex(v1->point()));
					glVertex(to_geex(v2->point()));
				}
			}
		}
		glEnd();
		std::cout << count/2 << " short edges." << std::endl;
	}

	void DelaunayGraphics::draw_isolate_long_edge()
	{
		std::set<Delaunay::Vertex_handle> incident_vertices; //vertices incident to a long edge
		int count = 0;
		glLineWidth(3) ;
		glColor3f(0.0, 0.0, 0.5) ;
		glBegin(GL_LINES) ;		
		FOR_EACH_FACE_DT(Delaunay, delaunay_, f) {
			if (delaunay_->is_infinite(f)) continue;
			for (int i=0; i<3; ++i)
			{
				Delaunay::Vertex_handle v1 = f->vertex((i+1)%3);
				Delaunay::Vertex_handle v2 = f->vertex((i+2)%3);
				if (delaunay_->is_long_edge(v1, v2))
				{
					if (incident_vertices.find(v1) != incident_vertices.end() ||
						incident_vertices.find(v2) != incident_vertices.end()) continue; //isolate edge condition
					count++;
					glVertex(to_geex(v1->point()));
					glVertex(to_geex(v2->point()));
					incident_vertices.insert(v1);
					incident_vertices.insert(v2);
				}
			}
		}
		glEnd();
		std::cout << count << " long edges." << std::endl;
	}

	void DelaunayGraphics::draw_isolate_short_edge()
	{
		std::set<Delaunay::Vertex_handle> incident_vertices; //vertices incident to a short edge
		int count = 0;
		glLineWidth(3) ;
		glColor3f(0.5, 0.0, 0.0) ;
		glBegin(GL_LINES) ;		
		FOR_EACH_FACE_DT(Delaunay, delaunay_, f) {
			if (delaunay_->is_infinite(f)) continue;
			for (int i=0; i<3; ++i)
			{
				Delaunay::Vertex_handle v1 = f->vertex((i+1)%3);
				Delaunay::Vertex_handle v2 = f->vertex((i+2)%3);
				if (delaunay_->is_short_edge(v1, v2))
				{
					if (incident_vertices.find(v1) != incident_vertices.end() ||
						incident_vertices.find(v2) != incident_vertices.end()) continue; //isolate edge condition
					count++;
					glVertex(to_geex(v1->point()));
					glVertex(to_geex(v2->point()));
					incident_vertices.insert(v1);
					incident_vertices.insert(v2);
				}
			}
		}
		glEnd();
		std::cout << count << " short edges." << std::endl;
	}

	void DelaunayGraphics::draw_separatrix() {
		delaunay_->calc_vertex_singularity();

		//collect irregular vertices && clear separatrix mark
		std::vector<Delaunay::Vertex_handle> irregular_vertices;
		FOR_EACH_VERTEX_DT(Delaunay, delaunay_, it) {
			it->separatrix_in_degree = 0;
			if (it->singularity != 0)
			{
				irregular_vertices.push_back(it);
			}
		}

		//collect separatrix edges
		typedef std::pair<Delaunay::Vertex_handle, Delaunay::Face_handle> dirEdge;  //<v,f> : Edge(v, incident_vertices(v,f))

		std::vector<dirEdge> separatrix0;
		for (int i=0; i<irregular_vertices.size(); ++i)
		{
			Delaunay::Vertex_handle v = irregular_vertices[i];
			int degree = delaunay_->dual_facet_degree(v);
			Delaunay::Face_circulator f = delaunay_->incident_faces(v);
			std::queue<dirEdge> Q;
			do 
			{
				Q.push(std::make_pair(v, f));
				f++;
			} while (f != delaunay_->incident_faces(v));

			while (!Q.empty())
			{
				dirEdge top = Q.front();
				Q.pop();
				separatrix0.push_back(top);

				Delaunay::Vertex_handle v1 = top.first;
				Delaunay::Face_handle   f1 = top.second;
				Delaunay::Vertex_handle v2 = delaunay_->incident_vertices(v1, f1);
				v2->separatrix_in_degree += 1;

				if (!v2->dual_intersects_boundary 
					&& v2->singularity == 0)
					//&& v2->separatrix_in_degree < 2)
				{
					Delaunay::Face_circulator fc = delaunay_->incident_faces(v2, f1);
					fc--; fc--;
					Q.push(std::make_pair(v2, fc));
				}
			}
		}

		//
		//std::vector<dirEdge> separatrix1;
		//for (int i=0; i<irregular_vertices.size(); ++i)
		//{
		//	Delaunay::Vertex_handle v = irregular_vertices[i];
		//	int degree = delaunay_->dual_facet_degree(v);
		//	Delaunay::Face_circulator f = delaunay_->incident_faces(v);
		//	std::queue<dirEdge> Q;
		//	do 
		//	{
		//		Q.push(std::make_pair(v, f));
		//		f++;
		//	} while (f != delaunay_->incident_faces(v));

		//	while (!Q.empty())
		//	{
		//		dirEdge top = Q.front();
		//		Q.pop();
		//		separatrix1.push_back(top);

		//		Delaunay::Vertex_handle v1 = top.first;
		//		Delaunay::Face_handle   f1 = top.second;
		//		Delaunay::Vertex_handle v2 = delaunay_->incident_vertices(v1, f1);
		//		//v2->separatrix_in_degree += 1;

		//		if (!v2->dual_intersects_boundary 
		//			&& v2->singularity == 0
		//			&& v2->separatrix_in_degree < 2)
		//		{
		//			Delaunay::Face_circulator fc = delaunay_->incident_faces(v2, f1);
		//			fc--; fc--;
		//			Q.push(std::make_pair(v2, fc));
		//		}
		//	}
		//}

		//draw separatrix
		glLineWidth(3);
		glBegin(GL_LINES);
		glColor3f(1, 0, 0);
		for(int i=0; i<separatrix0.size(); ++i) {
			dirEdge& edge = separatrix0[i];
			Delaunay::Vertex_handle v1 = edge.first;
			Delaunay::Vertex_handle v2 = delaunay_->incident_vertices(v1, edge.second);
			glVertex2f(v1->point().x(), v1->point().y());
			glVertex2f(v2->point().x(), v2->point().y());
		}
		//
		/*glColor3f(0, 0, 1);
		for(int i=0; i<separatrix1.size(); ++i) {
		dirEdge& edge = separatrix1[i];
		Delaunay::Vertex_handle v1 = edge.first;
		Delaunay::Vertex_handle v2 = delaunay_->incident_vertices(v1, edge.second);
		glVertex2f(v1->point().x(), v1->point().y());
		glVertex2f(v2->point().x(), v2->point().y());
		}*/
		glEnd();

		//
		std::cout << "separatrix length = " << separatrix0.size() << std::endl;

	}


	void DelaunayGraphics::draw_separatrix_period() {
		delaunay_->calc_separatrix_graph();
		const std::vector<Separatrix>& S = delaunay_->separatrix_graph_;

		glLineWidth(3);
		glBegin(GL_LINES);
		glColor3f(1, 0, 0);
		for (int i=0; i<S.size(); ++i)
		{
			draw_Separatrix(S[i]);
		}
		glEnd();
	}

	void DelaunayGraphics::draw_edge_stretch() {
		delaunay_->calc_edge_stretch();
		delaunay_->calc_vertex_singularity();
		delaunay_->calc_face_group_mark();
		int hmark = delaunay_->face_group_size();
		std::cout << "hmark = " << hmark << std::endl;

		//calc the most stretched edge for each face group
		std::vector<EdgeValue> isolate_edges(hmark, std::make_pair(triEdge(), 0));
		FOR_EACH_FACE_DT(Delaunay, delaunay_, it) {
			int gmark = it->group_mark;
			if (gmark != -1)
			{
				for (int i=0; i<3; ++i)
				{
					Delaunay::Vertex_handle v1 = it->vertex((i+1)%3);
					Delaunay::Vertex_handle v2 = it->vertex((i+2)%3);

					//
					if (v1->singularity != 0 
						|| (v1->dual_intersects_boundary && v2->dual_intersects_boundary)) {

							double stretch = it->edge_stretch[i]; 

							if (!isolate_edges[gmark].first.is_valid()
								|| (abs(stretch) > abs(isolate_edges[gmark].second)))
							{
								isolate_edges[gmark] = std::make_pair(triEdge(it, i), stretch);
							} 
					}
				}
			}
		}

		//draw grouped faces
		glPolygonMode(GL_FRONT_AND_BACK ,GL_FILL);
		glBegin(GL_TRIANGLES);
		FOR_EACH_FACE_DT(Delaunay, delaunay_, it) {
			int gmark = it->group_mark;
			if (gmark != -1 && gmark <= stretch_all_rank_)
			{

				glColor3f(gmark * 1.0/hmark, 1 - gmark*1.0/hmark, 0);
				for(int i=0; i<3; ++i) {
					Delaunay::Vertex_handle v = it->vertex(i);
					vec2 p = to_geex(v->point());
					glVertex2f(p.x, p.y);
				}
			}
		}
		glEnd();

		//collect valid stretched edges
		std::priority_queue<EdgeValue, std::vector<EdgeValue>, EdgeValue_greater> short_edges;
		std::priority_queue<EdgeValue, std::vector<EdgeValue>, EdgeValue_less> long_edges;
		FOR_EACH_FACE_DT(Delaunay, delaunay_, it) {
			if (delaunay_->is_infinite(it)) continue;
			if (delaunay_->has_irregular(it)) {
				for (int i=0; i<3; ++i)
				{
					Delaunay::Vertex_handle v1 = it->vertex((i+1)%3);
					Delaunay::Vertex_handle v2 = it->vertex((i+2)%3);

					//
					if (v1->singularity != 0 
						|| (v1->dual_intersects_boundary && v2->dual_intersects_boundary)) {

							double stretch = it->edge_stretch[i]; 

							if (stretch > 0)
							{
								long_edges.push(std::make_pair(triEdge(it, i), stretch));
							}
							if (stretch < 0)
							{
								short_edges.push(std::make_pair(triEdge(it, i), stretch));
							}
					}
				}
			}
		}

		//Covert to sorted vector of unique edges
		///long edges
		std::vector<EdgeValue> longVec;
		if (!long_edges.empty()) {
			longVec.push_back(long_edges.top());
			long_edges.pop();
		}
		while (!long_edges.empty())
		{
			if (long_edges.top().first != longVec.back().first)
			{
				longVec.push_back(long_edges.top());
			}
			long_edges.pop();
		}
		//short edges
		std::vector<EdgeValue> shortVec;
		if(!short_edges.empty()) {
			shortVec.push_back(short_edges.top());
			short_edges.pop();
		}
		while (!short_edges.empty())
		{
			if (short_edges.top().first != shortVec.back().first)
			{
				shortVec.push_back(short_edges.top());
			}
			short_edges.pop();
		}

		//Draw Long Edges
		if (!longVec.empty() && stretch_long_rank_ >= 0)
		{
			//calc max stretch
			double max_stretch =  longVec[0].second;
			if (stretch_long_rank_ < longVec.size())
			{
				max_stretch = longVec[stretch_long_rank_].second;
				std::cout << "stretch: max[" << stretch_long_rank_ << "] " << max_stretch;
			}
			//draw
			glLineWidth(2.5) ;
			glBegin(GL_LINES);
			//draw front edges
			int front_size = (1+stretch_long_rank_ < longVec.size() ? (1+stretch_long_rank_) : longVec.size());
			for (int i=0; i<front_size; ++i)
			{
				glColor3f(0, 1, 1);
				draw_triEdge(longVec[i].first);
			}
			//draw rest edges
			for (int i=front_size; i<longVec.size(); ++i)
			{
				EdgeValue& ev = longVec[i];
				triEdge& e = ev.first;
				double stretch = ev.second;
				glColor3f(0, 0, stretch/max_stretch);
				draw_triEdge(e);
			}
			glEnd();

		}

		//Draw Short Edges
		if (!longVec.empty() && stretch_short_rank_ >= 0)
		{
			//calc min stretch
			double min_stretch =  shortVec[0].second;
			if (stretch_short_rank_ < shortVec.size())
			{
				min_stretch = shortVec[stretch_short_rank_].second;
				std::cout << " min[" << stretch_short_rank_ << "] " << min_stretch << std::endl;
			}
			//draw
			glLineWidth(2.5) ;
			glBegin(GL_LINES);
			//draw front edges
			int front_size = (1+stretch_short_rank_ < shortVec.size() ? (1+stretch_short_rank_) : shortVec.size());
			for (int i=0; i<front_size; ++i)
			{
				glColor3f(1, 0, 1);
				draw_triEdge(shortVec[i].first);
			}
			//draw rest edges
			for (int i=front_size; i<shortVec.size(); ++i)
			{
				EdgeValue& ev = shortVec[i];
				triEdge& e = ev.first;
				double stretch = ev.second;
				glColor3f(stretch/min_stretch, 0, 0);
				draw_triEdge(e);
			}
			glEnd();
		}

		//All Edges
		//merge short and long edges
		// 		std::vector<EdgeValue> allVec;
		// 		std::vector<EdgeValue>::iterator lIt = longVec.begin(), sIt = shortVec.begin();
		// 		while(true) {
		// 			if (lIt == longVec.end())
		// 			{
		// 				while (sIt != shortVec.end())
		// 				{
		// 					allVec.push_back(*sIt);
		// 					++sIt;
		// 				}
		// 				break;
		// 			}
		// 			if (sIt == shortVec.end())
		// 			{
		// 				while (lIt != longVec.end())
		// 				{
		// 					allVec.push_back(*lIt);
		// 					++lIt;
		// 				}
		// 				break;
		// 			}
		// 
		// 			if (lIt->second > abs(sIt->second))
		// 			{
		// 				allVec.push_back(*lIt);
		// 				++lIt;
		// 			}
		// 			else {
		// 				allVec.push_back(*sIt);
		// 				++sIt;
		// 			}
		// 		}
		/*if (!allVec.empty() && stretch_all_rank_ >= 0)
		{
		int front_size = (1+stretch_all_rank_ < allVec.size()) ? (1+stretch_all_rank_) : allVec.size();
		std::cout << "stretch all[" << front_size -1 << "] " << allVec[front_size-1].second << std::endl;
		glLineWidth(2.5);
		glBegin(GL_LINES);
		for(int i=0; i<front_size; ++i) {
		EdgeValue& ev = allVec[i];
		if (ev.second > 0)
		{
		glColor3f(0, 1, 0.8);
		}
		else {
		glColor3f(1, 0.5, 0);
		}
		draw_triEdge(ev.first);
		}
		for(int i=front_size; i<allVec.size(); ++i) {
		glColor3f(0, 0, 0);
		draw_triEdge(allVec[i].first);
		}
		glEnd();
		}*/

		//draw isolated stretched edges
		if (!isolate_edges.empty() || stretch_all_rank_ >= 0)
		{
			int front_size = (1+stretch_all_rank_ < isolate_edges.size()) ? (1+stretch_all_rank_) : isolate_edges.size();
			glLineWidth(2.5);
			glBegin(GL_LINES);
			for (int i=0; i<front_size; ++i)
			{
				triEdge& e = isolate_edges[i].first;
				double stretch = isolate_edges[i].second;

				if (stretch > 0)
				{
					glColor3f(0, 1, 0.8);
				} 
				else
				{
					glColor3f(1, 0.5, 0);
				}
				draw_triEdge(e);
			}
			glEnd();
		}

	}

	void DelaunayGraphics::draw_edge_stretch_period() {
		delaunay_->calc_edge_stretch_period();
		delaunay_->calc_vertex_singularity_period();
		delaunay_->calc_face_group_mark_period();
		int hmark = delaunay_->face_group_size();
		std::cout << "hmark = " << hmark << std::endl;

		//calc the most stretched edge for each face group
		std::vector<EdgeValue> isolate_edges(hmark, std::make_pair(triEdge(), 0));
		FOR_EACH_FACE_DT(Delaunay, delaunay_, it) {
			int gmark = it->group_mark;
			if (gmark != -1)
			{
				for (int i=0; i<3; ++i)
				{
					Delaunay::Vertex_handle v1 = it->vertex((i+1)%3);
					Delaunay::Vertex_handle v2 = it->vertex((i+2)%3);

					//
					if (v1->singularity != 0 
						/*|| (v1->dual_intersects_boundary && v2->dual_intersects_boundary)*/) {

							double stretch = it->edge_stretch[i]; 
							//test: topo bias
							// 							int d1 = delaunay_->dual_facet_degree_period(v1);
							// 							int d2 = delaunay_->dual_facet_degree_period(v2);
							// 							int topo_bias = abs(d1 + d2 - 12);
							// 							stretch *= topo_bias;
							//

							if (!isolate_edges[gmark].first.is_valid()
								|| (abs(stretch) > abs(isolate_edges[gmark].second)))
							{
								isolate_edges[gmark] = std::make_pair(triEdge(it, i), stretch);
							} 
					}
				}
			}
		}

		//draw grouped faces
		glPolygonMode(GL_FRONT_AND_BACK ,GL_FILL);
		glBegin(GL_TRIANGLES);
		FOR_EACH_FACE_DT(Delaunay, delaunay_, it) {
			int gmark = it->group_mark;
			if (gmark != -1 && gmark <= stretch_all_rank_)
			{

				glColor3f(gmark * 1.0/hmark, 1 - gmark*1.0/hmark, 0);
				for(int i=0; i<3; ++i) {
					Delaunay::Vertex_handle v = it->vertex(i);
					vec2 p = to_geex(v->point());
					glVertex2f(p.x, p.y);
				}
			}
		}
		glEnd();

		//collect valid stretched edges
		std::priority_queue<EdgeValue, std::vector<EdgeValue>, EdgeValue_greater> short_edges;
		std::priority_queue<EdgeValue, std::vector<EdgeValue>, EdgeValue_less> long_edges;
		FOR_EACH_FACE_DT(Delaunay, delaunay_, it) {
			if (delaunay_->is_infinite(it)) continue;
			if (delaunay_->has_irregular(it)) {
				for (int i=0; i<3; ++i)
				{
					Delaunay::Vertex_handle v1 = it->vertex((i+1)%3);
					Delaunay::Vertex_handle v2 = it->vertex((i+2)%3);

					//
					if (v1->singularity != 0 
						/*|| (v1->dual_intersects_boundary && v2->dual_intersects_boundary)*/) {

							double stretch = it->edge_stretch[i]; 

							if (stretch > 0)
							{
								long_edges.push(std::make_pair(triEdge(it, i), stretch));
							}
							if (stretch < 0)
							{
								short_edges.push(std::make_pair(triEdge(it, i), stretch));
							}
					}
				}
			}
		}

		//Covert to sorted vector of unique edges
		///long edges
		std::vector<EdgeValue> longVec;
		if (!long_edges.empty()) {
			longVec.push_back(long_edges.top());
			long_edges.pop();
		}
		while (!long_edges.empty())
		{
			if (long_edges.top().first != longVec.back().first)
			{
				longVec.push_back(long_edges.top());
			}
			long_edges.pop();
		}
		//short edges
		std::vector<EdgeValue> shortVec;
		if(!short_edges.empty()) {
			shortVec.push_back(short_edges.top());
			short_edges.pop();
		}
		while (!short_edges.empty())
		{
			if (short_edges.top().first != shortVec.back().first)
			{
				shortVec.push_back(short_edges.top());
			}
			short_edges.pop();
		}

		//Draw Long Edges
		if (!longVec.empty() && stretch_long_rank_ >= 0)
		{
			//calc max stretch
			double max_stretch =  longVec[0].second;
			if (stretch_long_rank_ < longVec.size())
			{
				max_stretch = longVec[stretch_long_rank_].second;
				std::cout << "stretch: max[" << stretch_long_rank_ << "] " << max_stretch;
			}
			//draw
			glLineWidth(2.5) ;
			glBegin(GL_LINES);
			//draw front edges
			int front_size = (1+stretch_long_rank_ < longVec.size() ? (1+stretch_long_rank_) : longVec.size());
			for (int i=0; i<front_size; ++i)
			{
				glColor3f(0, 1, 1);
				draw_triEdge(longVec[i].first);
			}
			//draw rest edges
			for (int i=front_size; i<longVec.size(); ++i)
			{
				EdgeValue& ev = longVec[i];
				triEdge& e = ev.first;
				double stretch = ev.second;
				glColor3f(0, 0, stretch/max_stretch);
				draw_triEdge(e);
			}
			glEnd();

		}

		//Draw Short Edges
		if (!longVec.empty() && stretch_short_rank_ >= 0)
		{
			//calc min stretch
			double min_stretch =  shortVec[0].second;
			if (stretch_short_rank_ < shortVec.size())
			{
				min_stretch = shortVec[stretch_short_rank_].second;
				std::cout << " min[" << stretch_short_rank_ << "] " << min_stretch << std::endl;
			}
			//draw
			glLineWidth(2.5) ;
			glBegin(GL_LINES);
			//draw front edges
			int front_size = (1+stretch_short_rank_ < shortVec.size() ? (1+stretch_short_rank_) : shortVec.size());
			for (int i=0; i<front_size; ++i)
			{
				glColor3f(1, 0, 1);
				draw_triEdge(shortVec[i].first);
			}
			//draw rest edges
			for (int i=front_size; i<shortVec.size(); ++i)
			{
				EdgeValue& ev = shortVec[i];
				triEdge& e = ev.first;
				double stretch = ev.second;
				glColor3f(stretch/min_stretch, 0, 0);
				draw_triEdge(e);
			}
			glEnd();
		}

		//draw isolated stretched edges
		if (!isolate_edges.empty() || stretch_all_rank_ >= 0)
		{
			int front_size = (1+stretch_all_rank_ < isolate_edges.size()) ? (1+stretch_all_rank_) : isolate_edges.size();
			glLineWidth(2.5);
			glBegin(GL_LINES);
			for (int i=0; i<front_size; ++i)
			{
				triEdge& e = isolate_edges[i].first;
				double stretch = isolate_edges[i].second;

				if (stretch > 0)
				{
					glColor3f(0, 1, 0.8);
				} 
				else
				{
					glColor3f(1, 0.5, 0);
				}
				draw_triEdge(e);
			}
			glEnd();
		}

	}

	void DelaunayGraphics::draw_selected() {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glPointSize(10); ///test
		glLineWidth(3);

		FOR_EACH_FACE_DT(Delaunay, delaunay_, it) {
			if (delaunay_->is_infinite(it)) continue;

			//Vertex
			std::vector<Point> points;
			glColor3f(0.5, 0.0, 1.0);
			glBegin(GL_POINTS);
			for (int i=0; i<3; ++i)
			{
				Delaunay::Vertex_handle v = it->vertex(i);
				points.push_back(v->point());
				if (v->selected)
				{
					glVertex2f(v->point().x(), v->point().y());
				}
			}
			glEnd();

			//Edge
			glColor3f(0.0, 0.0, 0.0);
			glBegin(GL_LINES);
			for (int i=0; i<3; ++i)
			{
				if (it->selected_edge[i])
				{
					int j1 = (i+1) % 3;
					int j2 = (i+2) % 3;
					glVertex2f(points[j1].x(), points[j1].y());
					glVertex2f(points[j2].x(), points[j2].y());
				}

			}
			glEnd();

			//Face
			if (it->selected)
			{
				glColor3f(1.0, 0.5, 0.0);
				glBegin(GL_TRIANGLES);
				for (int i=0; i<3; ++i)
				{
					glVertex2f(points[i].x(), points[i].y());
				}
				glEnd();
			}

		}

	}

	void DelaunayGraphics::draw_selected_separatrix() {
		if (delaunay_->period()) delaunay_->calc_vertex_singularity_period();
		else delaunay_->calc_vertex_singularity();

		int num = delaunay_->selected_separatrix_.size();
		glLineWidth(4);
		glBegin(GL_LINES);
		if (num == 1)
		{
			draw_Separatrix(*(delaunay_->selected_separatrix_.begin()), 160);
			std::cout << "separatrix length = " << delaunay_->selected_separatrix_.begin()->length() << std::endl;
		}
		else {
			int i=0;
			const vec3 c1(120, 1.0, 1.0);
			const vec3 c2(0.0, 1.0, 1.0);
			for (auto it = delaunay_->selected_separatrix_.begin(); it!= delaunay_->selected_separatrix_.end(); ++it)
			{
				double m = i * 1.0 / num;
				vec3 mix_hsv = mix(c1, c2, m);
				vec3 mix_rgb;
				hsv2rgb(mix_hsv, mix_rgb);
				glColor(mix_rgb);
				draw_Separatrix(*it);
				++i;
			}
		}
		glEnd();
	}

	void DelaunayGraphics::rotate_vector(vec2& v, double angle) {
		double c = cos(angle);
		double s = sin(angle);
		double new_x = c * v.x - s * v.y;
		double new_y = s * v.x + c * v.y;
		v = vec2(new_x, new_y);
	}

	void DelaunayGraphics::draw_vector(const vec2& v, const vec2& start) {
		double len = v.length();
		if (len == 0) return;
		float arror_size = 0.15;
		vec2 arror1 = v;
		rotate_vector(arror1, 0.9 * M_PI);
		arror1 *= arror_size;
		vec2 arror2 = v;
		rotate_vector(arror2, -0.9 * M_PI);
		arror2 *= arror_size;

		glBegin(GL_LINES) ;
		vec2 end = start + v;
		glVertex2f(start.x, start.y);
		glVertex2f(end.x, end.y);
		vec2 arror1_end = end + arror1;
		glVertex2f(end.x, end.y);
		glVertex2f(arror1_end.x, arror1_end.y);
		vec2 arror2_end = end + arror2;
		glVertex2f(end.x, end.y);
		glVertex2f(arror2_end.x, arror2_end.y);
		glEnd();
	}

	void DelaunayGraphics::draw_primal_triangle(Delaunay::Face_handle f) {
		glBegin(GL_TRIANGLES);
		for (int i=0; i<3; ++i)
		{
			glVertex2f(f->vertex(i)->point().x(), f->vertex(i)->point().y());
		}
		glEnd();
	}
	//triEdge
	void DelaunayGraphics::draw_triEdge(const triEdge& e) {
		glVertex(to_geex(e.start()->point()));
		glVertex(to_geex(e.end()->point()));
	}
	//dirEdge
	void DelaunayGraphics::draw_dirEdge(const dirEdge& e) {
		glVertex(to_geex(e.start()->point()));
		glVertex(to_geex(e.end()->point()));
	}
	//Separatrix
	void DelaunayGraphics::draw_Separatrix(const Separatrix& s) {
		traceEdge (*next_ptr)(const traceEdge&) = traceEdge::next;
		if (delaunay_->period()) next_ptr = traceEdge::next_period;

		traceEdge It = s.first_edge();
		int vstart = s.first_vertex()->index;
		//first edge
		if (!It.is_null())
		{
			draw_dirEdge(It);
			It = next_ptr(It);
		}
		//rest edges
		while (!It.is_null() && (It.start()->index != vstart))
		{
			draw_dirEdge(It);
			It = next_ptr(It);
		}
	}

	//draw separatrix with changing color
	void DelaunayGraphics::draw_Separatrix(const Separatrix& s, float color_h) {
		double color_v = 1.0;
		double dv = 1.0 / s.length();

		traceEdge (*next_ptr)(const traceEdge&) = traceEdge::next;
		if (delaunay_->period()) next_ptr = traceEdge::next_period;

		traceEdge It = s.first_edge();
		int vstart = s.first_vertex()->index;
		//first edge
		if (!It.is_null())
		{
			//gradually changing color
			vec3 hsv(color_h, 1, color_v);
			vec3 rgb;
			hsv2rgb(hsv, rgb);
			glColor(rgb);
			color_v -= dv;
			//
			draw_dirEdge(It);
			It = next_ptr(It);
		}
		//rest edges
		while (!It.is_null() && (It.start()->index != vstart))
		{
			//gradually changing color
			vec3 hsv(color_h, 1, color_v);
			vec3 rgb;
			hsv2rgb(hsv, rgb);
			glColor(rgb);
			color_v -= dv;
			//
			draw_dirEdge(It);
			It = next_ptr(It);
		}

	}
	///dxy add end
}
