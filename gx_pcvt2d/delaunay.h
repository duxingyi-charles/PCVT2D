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

#ifndef __DELAUNAY__
#define __DELAUNAY__

//#define NOMINMAX 

#include "delaunay_cgal.h"
#include "polygons.h"
#include <GL/gl.h>
///dxy
#include <queue>
#include "irreg_editor.h"
///

namespace Geex {

	///dxy move from delaunay_cvt
	enum PVDMode { NO_COPY, FULL_COPY, MIN_COPY } ;
#define PVDModeNames "NO_COPY | Full Copy" | "Min Copy"
	///

	enum SamplingMode { PD_CCVT, PD_VORONOI, PD_BOUNDARY, PD_DARTTHROW, PD_PENROSE } ;
#define SamplingModeNames "Vornoi | CCVT" | "Pure" | "DartThrow" | "Penrose"

	class DelaunayGraphics ;
	class DelaunayCVT ;
	class IrregEditor ;
	class Map ;
	class GridCell ;
	///dxy
	class triEdge;
	class dirEdge;
	class traceEdge;
	class Separatrix;
	class Pair57;
	enum PairMoveCode { None, Up_right, Up_left, Right, Left, Down_right, Down_left };  //test
	//
	class PeriodTriangulation;
	//______________________________________________________________________________________

	class SegmentDelaunay : public DelaunayBase {
		typedef DelaunayBase baseclass ;
	public:
		void insert(int id, const vec2& p1, const vec2& p2, double step_length) ;
		int locate(const vec2& p) ;
	protected:
		void insert(int id, const vec2& p) ;
	} ;

	class Delaunay : public DelaunayBase {
		typedef DelaunayBase baseclass ;

	public:
		Delaunay() ;
		~Delaunay() ;

		void save(const std::string& filename) ;
		void load(const std::string& filename) ;
		void load_boundary(const std::string& filename) ;
		void set_non_convex_mode(bool x) { non_convex_mode_ = x ; }

		void get_bbox(
			real& x_min, real& y_min, real& z_min,
			real& x_max, real& y_max, real& z_max
			) ;

		bool in_boundary(const vec2& p) { 
			return non_convex_mode_ ? boundary_.contains(p) : boundary_convex_.contains(p) ; 
		}

		GLboolean& insert_boundary() { return insert_boundary_ ; }
		//		GLboolean& insert_full_copy() { return insert_full_copy_ ; }
		//		GLboolean& period() { return period_ ; }

		void compute_inner_voronoi(Vertex_handle vh, std::vector<vec2>& poly) ;

		// ------------------------------------ Posson disk

		vec2 origin() { return 0.5*vec2(x_min_+y_min_, x_max_+y_max_) ; } 
		GLfloat& perturb() { return perturb_ ; }
		void do_perturb() ;
		void digree_hist(double& plt5, double& p5, double& p6, double& p7, double& pgt7) ;
		void edge_hist(double& p66, double& p56, double& p76, double& p55, double& p77, double& p57) ;
		std::vector<double>& edge_hist() { return edge_hist_ ; }


		void generate_poisson_disk() ;
		void sample_poisson_grid(double radius, int nmaxfail) ;
		void init_sample_grid(int res, double radius, int nmaxfail) ;
		void insert_grid_samples_period(int res) ;
		void insert_grid_samples_clipped(int res) ;
		//		bool sample_grid_gaps(int res, double radius) ;
		inline bool is_hit_period(vec2& p, int u, int v, int res, double sqr_radius) ;

		void sample_poisson(double radius=0.001, bool isTiled=false, bool maximize=1, int minMaxThrows=1000, 
			int multiplier=false, int relax=0, const char* method="Penrose") ;
		void set_vertices(std::vector<vec2>& points) ;
		double& sample_radius() { return sample_radius_ ; }
		SamplingMode& sampling_mode() { return sampling_mode_ ; }

		void fast_poisson_disk(double exclude_radius) ;
		Vertex_handle insert_vertex_periodic(vec2& p, double exclude_radius) ;
		void compute_voronoi_area(Delaunay::Vertex_handle& v, double exclude_radius) ;
		vec2 new_vertex_in_voronoi(Delaunay::Vertex_handle& v, double exclude_radius);
		void compute_gap_area(Delaunay::Vertex_handle& v, double exclude_radius) ;
		vec2 new_vertex_in_gap(Delaunay::Vertex_handle& v, double exclude_radius) ;
		bool on_convex_hull(Delaunay::Vertex_handle& v) ;
		void compute_triangle_voids() ;
		bool is_dual_gap(Delaunay::Face_handle f, double exclude_radius) ;
		void compute_void_region(Delaunay::Face_handle f, double exclude_radius, std::vector<vec2>& poly) ;
		void update_gaps(double exclude_radius) ;
		void compute_void_skeleton(Face_handle f, double exclude_radius) ;
		void compute_void_region(Face_handle f, double exclude_radius) ;
		std::vector<std::vector<vec2> >& void_skeleton() { return void_skeleton_ ; } 
		std::vector<std::vector<vec2> >& void_regions() { return void_regions_ ; } 
		void insert_point_poisson(vec2& p, double radius) ;
		void insert_copy_poisson(vec2& p, int index, double radius) ;
		void update_neighbor_cells(Vertex_handle& vh, double radius) ;

		bool pick_vertex(vec2& pt) ;
		bool pick_edge(vec2& pt) ;

		std::vector<Vertex_handle>&  vertices_sel() { return vertices_sel_ ; }// selected vertices 
		std::vector<Edge>& edges_sel() { return edges_sel_ ; }     // selected edges
		void clear_selection() {
			vertices_sel_.clear() ; 
			edges_sel_.clear() ; 
		}
		Map* map()             { return map_ ;    }
		IrregEditor* editor()  { return editor_ ; }
		bool& map_edit_mode()  { return map_edit_mode_ ; }
		// ------------------------------------ Delaunay 

		int nb_vertices() const { return all_vertices_.size() ; }
		void clear() ;
		void begin_insert() ;
		///dxy add
		Vertex_handle move_nearest_vertex_to(const vec2& p);
		///
		Vertex_handle insert(const vec2& p) ;
		Vertex_handle insert(const vec2& p, int index) ;
		Vertex_handle nearest(const vec2& p) ;
		void remove(const vec2& p) ;
		void end_insert(bool redraw = true) ;

		void insert_random_vertex() ;
		void insert_random_vertices(int nb) ;
		void insert_grid(int nb) ;

		std::vector<Vertex_handle>& vertices() { return all_vertices_ ; }

		// ------------------------------------ Periodic Delaunay
		void compute_pvd() ;
		void compute_pvd2() ;
		void clear_copies(bool redraw=false) ;
		void insert_copies(bool full_copy=false, bool redraw=false) ;
		void begin_insert_copies() ;
		void end_insert_copies(bool redraw=false) ;
		void insert_copies_full(bool redraw=false) ;
		void insert_copies_ring(bool redraw=false) ;
		void insert_copies_poisson(double radius, bool redraw=false) ;

		bool is_primary(Vertex_handle v) { return v->domain == 0 ; }
		bool neighbor_to_primary(Vertex_handle v) ;
		bool is_full_hex() ; // check pvd configuration
		int  nb_primary() ;
		vec2 copy_point(vec2& p, PolygonEdge& e) ; // 
		vec2 copy_point(vec2& p, PolygonVertex& e) ; // 
		int  domain_idx(vec2& p) ;
		vec2 mirror_vertex_point(vec2& p, int vidx) ; // 
		vec2 mirror_edge_point(vec2& p, int eidx) ; // 
		vec2 translate(int domain_idx, vec2 p) ; // translate periodic coordinate to Euclidean
		std::set<Vertex_handle>& mirrors(int pid) {
			return mirrors_[pid] ;
		}
		void compute_edge_length() ; // test whether the edges of a hex with same length
		std::map<int, std::set<Vertex_handle> >& mirrors() { return mirrors_ ; }

		//void save_vertices(const std::string& filename) ;
		//void load_vertices(const std::string& filename) ;

		// ----------------------------------- RVD 2D
		//		typedef std::vector<PolygonEdge> PolyLine ;
		std::vector<Polygon2> rvd_ ;
		void compute_rvd(bool closed=true) ;
		bool is_boundary_cell(Delaunay::Vertex_handle v) { return rvd_[v->index].size()>0 ; }

		// ------------------------------------ Delaunay combinatorics ----------------

		vec2 dual(Face_handle c) { return c->dual ; }
		bool dual_cell_intersects_boundary(Vertex_handle v) {
			return v->dual_intersects_boundary ;
		}
		bool dual_cell_infinite(Vertex_handle v) {
			return v->dual_infinite ;
		}

		Line<real> get_dual_line(Edge e) {
			if(dimension()==1)
				return median_line(to_geex(e.first->vertex(1)->point()), to_geex(e.first->vertex(0)->point())) ;
			else
				return median_line(
				to_geex(e.first->vertex(ccw(e.second))->point()),
				to_geex(e.first->vertex( cw(e.second))->point())
				) ;
		}

		Polygon2* dual_convex_clip(Vertex_handle v, bool close = true) ;
		Polygon2* dual_convex_clip(Polygon2& from, Vertex_handle v, bool close = true) ;
		int dual_facet_degree(Vertex_handle v, bool period=false) ;
		int dual_facet_degree_period(Vertex_handle v) ;

		// --------------

		SegmentDelaunay& segments() { return segments_ ; }

		///dxy add
		//predicate
		bool is_vertex57(const Vertex_handle& v1, const Vertex_handle& v2);
		//short long topo edit
		void split_long_edge();
		void collapse_short_edge();
		void split_isolate_long_edge();
		void collapse_isolate_short_edge();
		//
		void split_isolate_long_edge_period();
		void collapse_isolate_short_edge_period();
		//select
		enum SelectMode {VERTEX, EDGE, FACE, SEPARATRIX_SEGMENT };
		GLenum& select_mode() { return select_mode_ ; }
		void clear_select();
		void update_select(const vec2&);
		//selected edge edit
		void collect_selected_edge(std::vector<triEdge>&);
		void collect_edge_middle(const std::vector<triEdge>&, std::vector<vec2>&);
		void collect_edge_ends(const std::vector<triEdge>&, std::vector<vec2>&);
		void split_selected_edge_period();
		void collapse_selected_edge_period();
		void collect_selected_Pair57(std::vector<Pair57>&);
		void collect_Pair57_operation(std::vector<Pair57> &P, std::vector<PairMoveCode> &M);
		//move Pair57
		void move_Pair57(std::vector<Pair57>&, std::vector<PairMoveCode>&);
		//
		void split_selected_edge();
		void collapse_selected_edge();
		void flip_selected_edge();
		//topo delaunay optimize
		void topo_edit_via_map(); //test for topo_delaunay_optimize
		void flip_selected_edge_via_map();  //test for topo_delaunay_optimize
		GLboolean& use_topo_delaunay_optimize() { return use_topo_delaunay_optimize_; }
		//stretch topo edit
		void calc_edge_stretch();
		void calc_vertex_singularity();
		void calc_face_group_mark();  // faces into non-adjacent groups
		void stretch_topo_optimize();
		int face_group_size() { 
			calc_face_group_mark();
			return face_group_size_; 
		}
		bool has_irregular(Delaunay::Face_handle);
		bool has_primary(Face_handle);
		//
		void calc_edge_stretch_period();
		void calc_vertex_singularity_period();
		void calc_face_group_mark_period();
		void calc_vertex_group_mark();
		//separatrix
		void calc_separatrix_graph();
		//
		void stretch_topo_optimize_period();
		void stretch_topo_optimize_period_balanced();
		//
		void update_area(bool info = true);
		//new add
		bool is_long_edge(const Vertex_handle& v1, const Vertex_handle& v2) 
		{
			if (!(is_primary(v1) || is_primary(v2))) return false;
			int d1 = dual_facet_degree_period(v1);
			int d2 = dual_facet_degree_period(v2);
			return (d1 + d2 >= 14); //long edge condition
		}
		bool is_short_edge(const Vertex_handle& v1, const Vertex_handle& v2)
		{
			if (!(is_primary(v1) || is_primary(v2))) return false;
			int d1 = dual_facet_degree_period(v1);
			int d2 = dual_facet_degree_period(v2);
			return (d1 + d2 <= 10); //short edge condition
		}
		bool is_long_edge(const triEdge& e);
		bool is_short_edge(const triEdge& e);
		//insert/remove
		void insert_array(const std::vector<vec2>& arr) {
			for (int i=0; i<arr.size(); ++i) insert(arr[i]);
		}
		void remove_array(const std::vector<vec2>& arr) {
			for (int i=0; i<arr.size(); ++i) remove(arr[i]);
		}
		///dxy add end
		///dxy move from delaunay_cvt
		GLenum& pvd_mode() { return pvd_mode_ ; }
		GLboolean period() { return pvd_mode_ != NO_COPY ; }
		void get_primary_position(vec2& g);
		///
		//////////////////////////////////////////////////////////////////////////dxy test func
		void print_primary_count();
		void print_grad_count(int);
		void print_topo_info(int num);
		//////////////////////////////////////////////////////////////////////////

	protected:
		std::vector<Vertex_handle> all_vertices_ ;
		std::map<int, std::vector<Vertex_handle> > mirror_vertices_ ;
		bool non_convex_mode_ ;

		Polygon2 boundary_ ;
		Convex boundary_convex_ ;
		bool cached_bbox_ ;
		double x_min_, y_min_, x_max_, y_max_ ;

		Polygon2 ping_ ;
		Polygon2 pong_ ;

		SegmentDelaunay segments_ ;

		friend class DelaunayGraphics ;
		friend class DelaunayCVT ;

		bool opened_ ;
		GLboolean insert_boundary_ ;
		//	GLboolean period_ ;
		//		GLboolean insert_full_copy_ ;
		unsigned int nb_master_ ;
		const static vec2 v_offset_[4] ;
		const static vec2 e_offset_[4] ;
		const static vec2 domain_offset_[9] ;

		// used for poisson disk sampling
		GLfloat perturb_ ;
		std::vector<double> edge_hist_ ;
		int bin_size_ ;
		double sample_radius_ ; 
		GridCell **grid_ ;
		std::map<int, std::vector<Vertex_handle> > vertices_ ; 
		std::map<int, std::set<Vertex_handle> > mirrors_ ;
		SamplingMode sampling_mode_ ;
		std::vector<std::vector<vec2> > void_skeleton_ ;
		std::vector<std::vector<vec2> > void_regions_ ;

		// user interaction
		std::vector<Vertex_handle>  vertices_sel_ ; // selected vertices 
		std::vector<Delaunay::Edge> edges_sel_ ;     // selected edges

		bool        map_edit_mode_ ;
		IrregEditor *editor_ ;
		Map         *map_ ;

		///dxy add
		PeriodTriangulation *periodic_tri_;
		//
		bool is_edge_stretch_dirty_;
		bool is_vertex_singularity_dirty_;
		bool is_face_group_mark_dirty_;
		bool is_vertex_group_mark_dirty_;
		int face_group_size_;

		float total_area_; //used to normalize cvt and direction energy
		float average_area_;

		std::set<Separatrix> selected_separatrix_;

		bool is_separatrix_graph_dirty_;
		std::vector<Separatrix> separatrix_graph_;
		GLenum select_mode_;
		//
		GLboolean use_topo_delaunay_optimize_;
		///dxy move from delaunay_cvt
		GLenum pvd_mode_ ;
		///
		//dxy add: friends
		friend class dirEdge;
		friend class traceEdge;


	} ;

	//______________________________________________________________________________________

#define FOR_EACH_VERTEX_DT(TCLASS, INSTANCE, IT)                       \
	for(                                                            \
	TCLASS::Vertex_iterator IT = (INSTANCE)->vertices_begin() ; \
	IT != (INSTANCE)->vertices_end(); IT++                      \
	)

	///dxy add
#define FOR_EACH_PRIMARY_VERTEX_DT(TCLASS, INSTANCE, IT)                       \
	for(                                                            \
	TCLASS::Vertex_iterator IT = (INSTANCE)->vertices_begin() ; \
	IT != (INSTANCE)->vertices_end(); IT++                      \
	) if (INSTANCE->is_primary(IT))
	///

#define FOR_EACH_FACE_DT(TCLASS, INSTANCE, IT)                      \
	for(                                                         \
	TCLASS::Face_iterator IT = (INSTANCE)->faces_begin() ;   \
	IT != (INSTANCE)->faces_end(); IT++                      \
	)

#define FOR_EACH_FINITE_FACE_DT(TCLASS, INSTANCE, IT)                              \
	for(                                                                        \
	TCLASS::Finite_faces_iterator IT = (INSTANCE)->finite_faces_begin() ;   \
	IT != (INSTANCE)->finite_faces_end(); IT++                              \
	)

#define FOR_EACH_EDGE_DT(TCLASS, INSTANCE, IT)                      \
	for(                                                         \
	TCLASS::Edge_iterator IT = (INSTANCE)->edges_begin() ;   \
	IT != (INSTANCE)->edges_end(); IT++                      \
	)


	///dxy add: triEdge 
	class triEdge {
	public:
		triEdge(Delaunay::Face_handle f = nullptr, unsigned idx = 0) : face_(f), edge_index_(idx)
		{
			start_ = 0, end_ = 0;
			if (f!=nullptr && idx < 3)
			{
				start_ = f->vertex((idx + 1) % 3);
				end_ = f->vertex((idx + 2) % 3);
			}
		}
		Delaunay::Vertex_handle start() const { return start_; }
		Delaunay::Vertex_handle end() const { return end_; }
		Delaunay::Face_handle face() const { return face_; }
		vec2 middle_point() const { return 0.5 * (to_geex(start_->point()) + to_geex(end_->point())); }
		unsigned edge_index() const { return edge_index_; }
		bool operator==(const triEdge &rhs) const {
			return (this->start() == rhs.start() && this->end() == rhs.end())
				|| (this->start() == rhs.end() && this->end() == rhs.start());
		}
		bool operator!=(const triEdge &rhs) const {	return !(*this == rhs);	}
		bool is_valid() const { return (face_ != nullptr) && (edge_index_ <= 3); }
	private:
		Delaunay::Face_handle face_;
		unsigned edge_index_;
		Delaunay::Vertex_handle start_;
		Delaunay::Vertex_handle end_;
	};


	typedef std::pair<triEdge, double> EdgeValue;
	struct EdgeValue_greater{  
		bool operator ()(EdgeValue &a, EdgeValue &b){  
			return a.second > b.second;  
		}  
	}; 
	struct EdgeValue_less{  
		bool operator ()(EdgeValue &a, EdgeValue &b){  
			return a.second < b.second;  
		}  
	};

	//dirEdge
	class dirEdge
	{
	public:
		dirEdge(Delaunay::Vertex_handle v = nullptr, Delaunay::Face_handle f = nullptr);
		static bool is_null(const dirEdge& e) { return (e.v_ == nullptr || e.f_ == nullptr); }
		static bool is_valid(const Delaunay::Vertex_handle& v, const Delaunay::Face_handle& f) {
			if (v == nullptr || f == nullptr) return false;
			for (int i=0; i<3; ++i) if (v == f->vertex(i)) return true;
			return false;
		}
		bool is_valid() const { return is_valid(v_, f_); }
		bool is_null() const { return is_null(*this)/*(v_ == nullptr || f_ == nullptr)*/; }

		bool operator==(const dirEdge& rhs) const { return (this->v_ == rhs.v_ && this->f_ == rhs.f_); }
		bool operator!=(const dirEdge& rhs) const { return !(*this == rhs); }

		double length() const {
			if (v_ != nullptr && end_ != nullptr)
			{
				vec2 p1 = to_geex(v_->point());
				vec2 p2 = to_geex(end_->point());
				return Geex::length(p1-p2);
			} 
			else return std::numeric_limits<double>::quiet_NaN();
		}

		vec2 get_vector() const {
			if (v_ != nullptr && end_ != nullptr)
			{
				vec2 p1 = to_geex(v_->point());
				vec2 p2 = to_geex(end_->point());
				return (p2 - p1);
			}
			else return vec2(0, 0);
		}

		static bool is_same_edge(const dirEdge& a, const dirEdge& b) {
			if (a.is_valid() && b.is_valid())
			{
				if ((a.start() == b.start() && a.end() == b.end())
					|| (a.start() == b.end() && a.end() == b.start())) return true;
			}
			return false;
		}

		Delaunay::Face_handle face() const { return f_; }
		Delaunay::Vertex_handle start() const { return v_; }
		Delaunay::Vertex_handle end() const { return end_; }
		vec2 middle_point() const { 
			vec2 mid = 0.5 * (to_geex(v_->point()) + to_geex(end_->point()));
			return mid;
		}

		static void set_delaunay(Delaunay* d) { delaunay_ = d; }

		static dirEdge reverse(const dirEdge& e) {
			Delaunay::Face_circulator fc = delaunay_->incident_faces(e.end_, e.f_);
			++fc;
			return dirEdge(e.end_, fc);
		}

		static dirEdge primary_edge(const dirEdge& e) {
			if (is_null(e)) return e;
			if (delaunay_->is_primary(e.v_) || delaunay_->is_primary(e.end_)) return e;
			else {
				vec2 direct = e.get_vector();
				Delaunay::Vertex_handle v1 = delaunay_->all_vertices_[e.v_->index];
				vec2 p1 = to_geex(v1->point());
				vec2 p2 = p1 + direct;
				Delaunay::Face_circulator fc = delaunay_->incident_faces(v1);
				double mindist = 100;
				Delaunay::Face_handle fmin;
				do 
				{
					Delaunay::Vertex_handle vc = delaunay_->incident_vertices(v1, fc);
					double dist = distance(p2, to_geex(vc->point()));
					if (dist < mindist)
					{
						mindist = dist;
						fmin = fc;
					}
					fc++;
				} while (fc != delaunay_->incident_faces(v1));
				return dirEdge(v1, fmin);
			}
		}

	protected:
		Delaunay::Vertex_handle v_;
		Delaunay::Face_handle f_;
		Delaunay::Vertex_handle end_;
		static Delaunay* delaunay_;
	};

	struct dirEdge_length_less
	{
		bool operator()(const dirEdge& a, const dirEdge& b) {
			return a.length() < b.length();
		}
	};

	//traceEdge
	class traceEdge : public dirEdge
	{
	public:
		traceEdge(Delaunay::Vertex_handle v = nullptr, Delaunay::Face_handle f = nullptr) : dirEdge(v, f) {};
		traceEdge(const dirEdge& e) : dirEdge(e) {};
		static traceEdge next(const traceEdge& e) {
			if (e.end_->singularity) // end vertex is singular, return null edge
			{
				return traceEdge(nullptr, nullptr);
			}
			Delaunay::Face_circulator fc = delaunay_->incident_faces(e.end_, e.f_);
			fc--; fc--;
			return traceEdge(e.end_, fc);
		}

		static traceEdge next_period(const traceEdge& e) {
			return primary_edge(next(e));
		}

		void single_trace(traceEdge& last, int& count) const;
		void double_trace(traceEdge& first, traceEdge& last, int& count) const;

	};



	inline void traceEdge::single_trace(traceEdge& last, int& count) const {
		traceEdge (*next_ptr)(const traceEdge&) = next;
		if (delaunay_->period()) next_ptr = next_period;

		count = 0;
		traceEdge It = *this;
		if (delaunay_->period()) It = primary_edge(*this);

		int idx_start = v_->index;
		do 
		{
			last = It;
			++count;
			It = next_ptr(It);
		} while (!is_null(It) && (It.v_->index != idx_start));
		last = reverse(last);
	}

	inline void traceEdge::double_trace(traceEdge& first, traceEdge& last, int& count) const {
		int cnt1;
		single_trace(last, cnt1);

		traceEdge e = reverse(*this);
		int cnt2;
		e.single_trace(first, cnt2);

		count = cnt1 + cnt2 - 1;
		if (last.v_ == this->v_ && last.v_ != first.v_) //no singularity on trace line
		{
			assert(cnt1 == cnt2); //otherwise, there is bug
			count = cnt1;
		}
	}


	//class Separatrix
	class Separatrix
	{
	public:
		Separatrix(const traceEdge& e1, const traceEdge& e2, int len) {
			ends_[0] = e1; ends_[1] = e2;
			length_ = len;
		}
		Separatrix(const traceEdge& e) {
			e.double_trace(ends_[0], ends_[1], length_);
		}

		static bool is_null(const Separatrix& s) { return s.ends_[0].is_null() || s.ends_[1].is_null(); }
		static bool is_valid(const Separatrix& s) {
			if (is_null(s)) return false;
			else
			{
				traceEdge last;
				int len;
				s.ends_[0].single_trace(last, len);
				if (last == s.ends_[1])
				{
					return true;
				}
				else return false;
			}
		}

		bool is_null() const { return is_null(*this); }
		bool is_valid() const { return is_valid(*this); }

		bool operator==(const Separatrix& rhs) const {
			return (ends_[0] == rhs.ends_[0] && ends_[1] == rhs.ends_[1])
				|| (ends_[0] == rhs.ends_[1] && ends_[1] == rhs.ends_[0]);
		}
		bool operator!=(const Separatrix& rhs) const { return !(*this == rhs); }

		traceEdge first_edge() const { return ends_[0]; }
		traceEdge last_edge() const { return ends_[1]; }
		int length() const { return length_; }

		Delaunay::Vertex_handle first_vertex() const { return ends_[0].start(); }
		Delaunay::Vertex_handle last_vertex() const { return ends_[1].start(); }

		//test: for set<less> of separatrix
		bool operator<(const Separatrix& rhs) const {
			if (this->length_ == rhs.length_)
			{
				double thislen = this->first_edge().length() + this->last_edge().length();
				double rhslen = rhs.first_edge().length() + rhs.last_edge().length();
				return thislen < rhslen;
			} 
			else
			{
				return this->length_ < rhs.length_;
			}
		}

	private:
		traceEdge ends_[2];
		int length_;
	};


	///class Pair57
	class Pair57 : public dirEdge {
	public:
		//enum Move_Code {None, Up_right, Up_left, Right, Left, Down_right, Down_left} ;

		Pair57(const Delaunay::Vertex_handle& v = nullptr, const Delaunay::Face_handle& f = nullptr); 
		Pair57(const dirEdge&);
		Pair57(const vec2& pos);

		static bool is_valid(const dirEdge& e) {
			return (e.is_valid() && delaunay_->is_vertex57(e.start(), e.end()));
		}
		static bool is_valid(const Delaunay::Vertex_handle& v, const Delaunay::Face_handle& f)
		{
			dirEdge e(v, f);
			return is_valid(e);
		}
		bool is_null() const { return (status_ == Null); }
		bool is_valid() const { return (status_ == Valid); }

		void set_position(const vec2& p);
		bool verify();

		bool virtual_move(PairMoveCode move, std::vector<vec2>& to_insert, std::vector<vec2>& to_remove, std::vector<vec2>& new_positions);


	private:
		static bool search(const vec2& pos, dirEdge& result)
		{
			if (Numeric::has_nan(pos)) return false; 
			else
			{
				bool found57 = false;
				Delaunay::Vertex_handle res_v = nullptr;
				Delaunay::Face_handle res_f = nullptr;

				vec2 p = pos;
				if (delaunay_->period()) delaunay_->get_primary_position(p);

				Delaunay::Face_handle f0 = delaunay_->locate(to_cgal(p));
				if (delaunay_->is_infinite(f0)) return false;
				else
				{ //search 1-ring of f0
					for (int i=0; i<3; ++i)
					{
						Delaunay::Vertex_handle v = f0->vertex(i);
						Delaunay::Face_circulator fc = delaunay_->incident_faces(v);
						do 
						{
							if (!delaunay_->is_infinite(fc))
							{
								for (int j=0; j<3; ++j)
								{
									Delaunay::Vertex_handle v1 = fc->vertex(j);
									Delaunay::Vertex_handle v2 = fc->vertex((j+1)%3);
									if (delaunay_->is_vertex57(v1, v2))
									{
										if (found57) {
											if (fc != res_f || v1 != res_v)
											{
												return false; //find more than one 5-7 pair!
											}
										}
										else
										{
											found57 = true;
											res_f = fc;
											res_v = v1;
										}
									}
								}
							}
							fc++;
						} while (fc != delaunay_->incident_faces(v));
					}
				}
				if (found57)
				{
					result = dirEdge(res_v, res_f);
					return true;
				}
				else return false;
			}
		}

	private:
		enum Status {Valid, Pos_Only, Null};
		int status_;
		vec2 position_;
	};

	//class Pair57
	inline Pair57::Pair57(const Delaunay::Vertex_handle& v, const Delaunay::Face_handle& f) 
		: dirEdge(v, f) {
			if (is_valid(v, f)) {
				position_ = middle_point();
				status_ = Valid;
			}
			else {
				status_ = Null;
			}
	}

	inline Pair57::Pair57(const dirEdge& e) : dirEdge(e) {
		if (is_valid(e))
		{
			position_ = middle_point();
			status_ = Valid;
		}
		else {
			status_ = Null;
		}
	}

	inline Pair57::Pair57(const vec2& pos) {
		if (Numeric::has_nan(pos))
		{
			status_ = Null;
		}
		else {
			position_ = pos;
			status_ = Pos_Only;
		}
	}

	inline void Pair57::set_position(const vec2& p) {
		if (Numeric::has_nan(p))
		{
			status_ = Null;
		}
		else {
			position_ = p;
			status_ = Pos_Only;
		}
	}

	inline bool Pair57::verify() {
		if (status_ == Pos_Only)
		{
			dirEdge e;
			if (search(position_, e))
			{
				*this = e;
				position_ = middle_point();
				status_ = Valid;
			}
			else {
				status_ = Null;
			}
		}
		return (status_ == Valid);
	}

	inline bool Pair57::virtual_move(PairMoveCode move, std::vector<vec2>& to_insert, std::vector<vec2>& to_remove, std::vector<vec2>& new_positions)
	{
		if (status_ == Null)
		{
			new_positions.push_back(vec2(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()));
			return false;
		}
		else if (status_ == Pos_Only)
		{
			verify();
			return virtual_move(move, to_insert, to_remove, new_positions);
		}
		else { //Valid

			/*Notice: 
			1. don't guarantee there is no duplicates in to_insert/to_remove
			2. don't deal with the situation where 2 or more 5-7 pairs are very near
			*/

			switch (move)
			{
			case None:
				new_positions.push_back(middle_point());  //Warning: what about period situation?
				break;
			case Up_right: 
				{
					vec2 old_p5 = to_geex(v_->point());
					Delaunay::Face_circulator fc = delaunay_->incident_faces(v_, f_);
					fc++;
					Delaunay::Vertex_handle v6 = delaunay_->incident_vertices(v_, fc);
					vec2 p6 = to_geex(v6->point());
					vec2 new_p7 = 0.5 * (old_p5 + p6);
					fc++;
					Delaunay::Vertex_handle new_v5 = delaunay_->incident_vertices(v_, fc);
					vec2 new_p5 = to_geex(new_v5->point());
					vec2 new_pos = 0.5 * (new_p5 + new_p7);
					if (delaunay_->period())
					{
						delaunay_->get_primary_position(new_p7);
						delaunay_->get_primary_position(old_p5);
						delaunay_->get_primary_position(p6);
						delaunay_->get_primary_position(new_pos);
					}
					to_insert.push_back(new_p7);
					to_remove.push_back(old_p5);
					to_remove.push_back(p6);
					new_positions.push_back(new_pos);
				}
				break;
			case Up_left:
				{
					vec2 old_p5 = to_geex(v_->point());
					Delaunay::Face_circulator fc = delaunay_->incident_faces(v_, f_);
					fc--;
					Delaunay::Vertex_handle v6 = delaunay_->incident_vertices(v_, fc);
					vec2 p6 = to_geex(v6->point());
					vec2 new_p7 = 0.5 * (old_p5 + p6);
					fc--;
					Delaunay::Vertex_handle new_v5 = delaunay_->incident_vertices(v_, fc);
					vec2 new_p5 = to_geex(new_v5->point());
					vec2 new_pos = 0.5 * (new_p5 + new_p7);
					if (delaunay_->period())
					{
						delaunay_->get_primary_position(new_p7);
						delaunay_->get_primary_position(old_p5);
						delaunay_->get_primary_position(p6);
						delaunay_->get_primary_position(new_pos);
					}
					to_insert.push_back(new_p7);
					to_remove.push_back(old_p5);
					to_remove.push_back(p6);
					new_positions.push_back(new_pos);
				}
				break;
			case Right:
				{
					vec2 old_p5 = to_geex(v_->point());
					Delaunay::Face_circulator fc = delaunay_->incident_faces(v_, f_);
					fc++;
					Delaunay::Vertex_handle v6_1 = delaunay_->incident_vertices(v_, fc);
					fc = delaunay_->incident_faces(v6_1, f_);
					fc++; fc++;
					Delaunay::Vertex_handle v6_2 = delaunay_->incident_vertices(v6_1, fc);
					Delaunay::Vertex_handle v6_3 = delaunay_->incident_vertices(v6_2, fc);
					vec2 p6_1 = to_geex(v6_1->point());
					vec2 p6_2 = to_geex(v6_2->point());
					vec2 p6_3 = to_geex(v6_3->point());
					vec2 new_p5 = 0.5 * (p6_1 + p6_3);
					vec2 new_pos = 0.5 * (new_p5 + p6_2);
					if (delaunay_->period())
					{
						delaunay_->get_primary_position(new_p5);
						delaunay_->get_primary_position(old_p5);
						delaunay_->get_primary_position(new_pos);
					}
					to_insert.push_back(new_p5);
					to_remove.push_back(old_p5);
					new_positions.push_back(new_pos);
				}
				break;
			case Left:
				{
					vec2 old_p5 = to_geex(v_->point());
					Delaunay::Face_circulator fc = delaunay_->incident_faces(v_, f_);
					fc--;
					Delaunay::Vertex_handle v6_1 = delaunay_->incident_vertices(v_, fc);
					fc = delaunay_->incident_faces(v6_1, fc);
					fc--; fc--;
					Delaunay::Vertex_handle v6_3 = delaunay_->incident_vertices(v6_1, fc);
					Delaunay::Vertex_handle v6_2 = delaunay_->incident_vertices(v6_3, fc);
					vec2 p6_1 = to_geex(v6_1->point());
					vec2 p6_2 = to_geex(v6_2->point());
					vec2 p6_3 = to_geex(v6_3->point());
					vec2 new_p5 = 0.5 * (p6_1 + p6_3);
					vec2 new_pos = 0.5 * (new_p5 + p6_2);
					if (delaunay_->period())
					{
						delaunay_->get_primary_position(new_p5);
						delaunay_->get_primary_position(old_p5);
						delaunay_->get_primary_position(new_pos);
					}
					to_insert.push_back(new_p5);
					to_remove.push_back(old_p5);
					new_positions.push_back(new_pos);
				}
				break;
			case Down_right:
				{
					Delaunay::Face_circulator fc = delaunay_->incident_faces(end_, f_);
					fc--;
					Delaunay::Vertex_handle v6 = delaunay_->incident_vertices(end_, fc);
					fc--;
					Delaunay::Vertex_handle new_v7 = delaunay_->incident_vertices(end_, fc);
					vec2 old_p7 = to_geex(end_->point());
					vec2 p6 = to_geex(v6->point());
					vec2 new_p7 = to_geex(new_v7->point());
					vec2 new_p5 = 0.5 * (old_p7 + p6);
					vec2 new_pos = 0.5 * (new_p5 + new_p7);
					if (delaunay_->period())
					{
						delaunay_->get_primary_position(new_p5);
						delaunay_->get_primary_position(new_pos);
					}
					to_insert.push_back(new_p5);
					new_positions.push_back(new_pos);
				}
				break;
			case Down_left:
				{
					Delaunay::Face_circulator fc = delaunay_->incident_faces(end_, f_);
					fc++; fc++; fc++;
					Delaunay::Vertex_handle v6 = delaunay_->incident_vertices(end_, fc);
					fc++;
					Delaunay::Vertex_handle new_v7 = delaunay_->incident_vertices(end_, fc);
					vec2 old_p7 = to_geex(end_->point());
					vec2 p6 = to_geex(v6->point());
					vec2 new_p7 = to_geex(new_v7->point());
					vec2 new_p5 = 0.5 * (old_p7 + p6);
					vec2 new_pos = 0.5 * (new_p5 + new_p7);
					if (delaunay_->period())
					{
						delaunay_->get_primary_position(new_p5);
						delaunay_->get_primary_position(new_pos);
					}
					to_insert.push_back(new_p5);
					new_positions.push_back(new_pos);
				}
				break;
			}
			return true;
		}
	}
	///
}

#endif
