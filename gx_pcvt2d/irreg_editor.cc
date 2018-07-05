#include "irreg_editor.h"
#include "delaunay.h"

#include <Geex/combinatorics/map.h>
#include <Geex/combinatorics/map_editor.h>
#include <Geex/combinatorics/map_builder.h>

namespace Geex {

	IrregEditor::IrregEditor(Delaunay *delaunay, Map *map) 
		: map_(map)
		, delaunay_(delaunay) {
			delaunay_to_map() ;
			set_target(map_) ;
			cur_f_ = nil ;
			edit_mode_ = EDGE_FLIP ;
	}

	void IrregEditor::delaunay_to_map() {
		std::cout << "delaunay_to_map, begin..." << std::endl;
		std::vector<Delaunay::Vertex_handle> vertices = delaunay_->vertices() ;
		std::map<int, std::set<Delaunay::Vertex_handle> >& mirrors = delaunay_->mirrors() ;
		std::map<Delaunay::Vertex_handle, int> indices ;
		
		std::cout << vertices.size() << " primary vertices." << std::endl;

		map_->clear() ;
		MapBuilder builder(map_) ;
		builder.begin_surface() ;

		int cur = 0 ;
		for(unsigned int i=0; i<vertices.size(); ++i) {
			vec2 v = to_geex(vertices[i]->point()) ;
			builder.add_vertex(vec3(v.x, v.y, 0)) ;
			indices[vertices[i]] = cur ;
			cur ++ ;

			std::set<Delaunay::Vertex_handle>& vmirror = mirrors[i] ;
			if(vmirror.size() > 0) {
				//std::cout << vmirror.size() << std::endl;  //vmirror.size()==8
				for(std::set<Delaunay::Vertex_handle>::iterator it=vmirror.begin(); it!=vmirror.end(); ++it) {
					vec2 v = to_geex((*it)->point()) ;
					builder.add_vertex(vec3(v.x, v.y, 0)) ;
					indices[*it] = cur ;
					cur ++ ;			
				}
			}
		}

		std::cout << map_->nb_vertices() << " map vertices." << std::endl;

		//dxy: vertex property
		if (delaunay_->pvd_mode() == FULL_COPY) {
// 			Property<Map::Vertex, bool> is_primary(map_->vertex_property_manager(), "is_primary", PropertyManager::CREATE);
// 			Property<Map::Vertex, int> primary_id(map_->vertex_property_manager(), "primary_id", PropertyManager::CREATE);  //CREATE cuase bug
			Property<Map::Vertex, bool> is_primary(map_->vertex_property_manager(), "is_primary", PropertyManager::FIND_OR_CREATE);
			Property<Map::Vertex, int> primary_id(map_->vertex_property_manager(), "primary_id", PropertyManager::FIND_OR_CREATE);
			std::cout << "vertex property created" << std::endl;

			for(unsigned int i=0; i<vertices.size(); ++i) {
				Map::Vertex *v = builder.vertex(indices[vertices[i]]);
				is_primary[v] = true;
				primary_id[v] = i;

				std::set<Delaunay::Vertex_handle>& vmirror = mirrors[i];
				if (vmirror.size() > 0) {
					for(std::set<Delaunay::Vertex_handle>::iterator it=vmirror.begin(); it!=vmirror.end(); ++it) {
						v = builder.vertex(indices[*it]);
						is_primary[v] = false;
						primary_id[v] = i;
					}
				}
			}

			std::cout <<"vertex property initialized" << std::endl;
		}
		//

		FOR_EACH_FINITE_FACE_DT(Delaunay, delaunay_, it) {
			builder.begin_facet() ;
			for(int i=0; i<3; ++i) {
				builder.add_vertex_to_facet(indices[it->vertex(i)]) ;
			}
			builder.end_facet() ;
		}

		builder.end_surface() ;	
	}

	void IrregEditor::map_to_delaunay() {
		if (delaunay_->pvd_mode() == FULL_COPY)
		{
			std::cout << "map to delaunay begins ... , " << map_->nb_vertices() << "map vertices" << std::endl;
			Property<Map::Vertex, bool> is_primary(map_->vertex_property_manager(), "is_primary", PropertyManager::FIND);
			delaunay_->clear() ;
			delaunay_->begin_insert() ;
			FOR_EACH_VERTEX(Map, map_, v) {
				if(is_primary[v]) {
					vec3 pt = v->point();
					vec2 p(pt.x, pt.y);
					delaunay_->get_primary_position(p);
					delaunay_->insert(p);
				}
			}
			delaunay_->end_insert(false) ;
			delaunay_->insert_copies(delaunay_->pvd_mode()==FULL_COPY, false) ;
			std::cout << delaunay_->vertices().size() << " primary vertices." << std::endl;
		}
	}

	void IrregEditor::clear_picked() {
		picked_verts_.clear() ;
		picked_edges_.clear() ;
		cur_f_ = nil ;
	}

	//dxy change
	void IrregEditor::do_editing() {
		switch(edit_mode_) {
		case EDGE_FLIP:
			for(auto it=picked_edges_.begin(); it!=picked_edges_.end(); ++it) {
				edge_flip(*it);
			}
			break ;
		case EDGE_COLLAPSE: {
			for(auto it=picked_edges_.begin(); it!=picked_edges_.end(); ++it) {
				Halfedge *h = *it ;
				Vertex *v = h->vertex() ;
				vec3 p = 0.5*(v->point()+h->opposite()->vertex()->point()) ;
				collapse_edge(*it) ;
				v->set_point(p) ;
			}
			break ;
							}
		case VERTEX_SPLIT: {
			std::vector<bool> visited(picked_edges_.size(), false);
			for(int i=0; i<picked_verts_.size(); ++i) {
				Vertex *v = picked_verts_[i];
				std::vector<Halfedge*> incident_edges;
				for(int j=0; j<picked_edges_.size(); ++j) {
					if(!visited[j]) {
						Halfedge *ej = picked_edges_[j];
						if (ej->vertex()==v) {
							incident_edges.push_back(ej->opposite());
							visited[j] = true;
						}
						else if (ej->opposite()->vertex()==v) {
							incident_edges.push_back(ej);
							visited[j] = true;
						}
					}
				}
				if (incident_edges.size()==2) {
					split_vertex(v, incident_edges[0], incident_edges[1]);
				}
			}
			break ;
						   }			
		case V4_SPLIT:
			//split_v4(picked_verts_[0], picked_edges_[0]) ;
			break ;
		case V4_CREATE:
			//create_v4(picked_edges_[0], picked_edges_[1]) ;
			break ;
		case VL7_SPLIT:
			break ;
		case VL7_CREATE:
			break ;
		default:
			std::cerr << "the editing mode is not implemented..." << std::endl ;
			break ;
		}

		clear_picked() ;		
	}

// 	void IrregEditor::do_editing() {
// 		switch(edit_mode_) {
// 		case EDGE_FLIP:
// 			edge_flip(picked_edges_[0]) ;
// 			break ;
// 		case EDGE_COLLAPSE: {
// 			Halfedge *h = picked_edges_[0] ;
// 			Vertex *v = h->vertex() ;
// 			vec3 p = 0.5*(v->point()+h->opposite()->vertex()->point()) ;
// 			collapse_edge(picked_edges_[0]) ;
// 			v->set_point(p) ;
// 			break ;
// 							}
// 		case VERTEX_SPLIT: {
// 			Vertex *v = picked_verts_[0] ;
// 			Halfedge *f1 = picked_edges_[0] ;
// 			Halfedge *g1 = picked_edges_[1] ;
// 			if(f1->vertex()==v) f1 = f1->opposite() ;
// 			if(g1->vertex()==v) g1 = g1->opposite() ;
// 			split_vertex(v, f1, g1) ;
// 			break ;
// 						   }			
// 		case V4_SPLIT:
// 			split_v4(picked_verts_[0], picked_edges_[0]) ;
// 			break ;
// 		case V4_CREATE:
// 			create_v4(picked_edges_[0], picked_edges_[1]) ;
// 			break ;
// 		case VL7_SPLIT:
// 			break ;
// 		case VL7_CREATE:
// 			break ;
// 		default:
// 			std::cerr << "the editing mode is not implemented..." << std::endl ;
// 			break ;
// 		}
// 
// 		clear_picked() ;		
// 	}

	//dxy change end

	void IrregEditor::split_v4(Vertex* v, Halfedge* h) {
		Halfedge* ih = h->vertex()==v ? h : h->opposite() ;
		Halfedge* flip = ih->prev() ;
		collapse_edge(ih->next()) ;
		edge_flip(flip->prev()) ;
	}

	void IrregEditor::create_v4(Halfedge* h1, Halfedge* h2) {
		gx_assert(h1->prev()->opposite()==h2->prev()) ;
		split_edge(h1->prev(), true) ;
	}

	bool IrregEditor::edge_flip(Map::Halfedge* h)	{
		if(!is_flippable(h))
			return false;

		Map::Halfedge* hopp = h->opposite();

		Map::Halfedge* h00 = h->prev();
		Map::Halfedge* h01 = h->next();

		Map::Halfedge* h10 = hopp->next();
		Map::Halfedge* h11 = hopp->prev();

		Map::Facet* f0 = h->facet();
		Map::Facet* f1 = hopp->facet();


		link(h, h11, 1);
		link(h11, h01, 1);
		link(h01, h, 1);

		link(hopp, h00, 1);
		link(h00, h10, 1);
		link(h10, hopp, 1);

		set_facet_on_orbit(h, f0);
		make_facet_key(h);


		set_facet_on_orbit(hopp, f1);
		make_facet_key(hopp);

		make_vertex_key(h, h10->vertex());
		make_vertex_key(hopp, h01->vertex());

		make_vertex_key(h00);
		make_vertex_key(h10);
		make_vertex_key(h01);
		make_vertex_key(h11);

		return true;
	}

	bool IrregEditor::is_flippable(Map::Halfedge* h)	{
		// the two edges involved in the flip

		vec3 plane_v0 = normalize(h->vertex()->point() - h->opposite()->vertex()->point());
		vec3 plane_v1 = normalize(h->opposite()->next()->vertex()->point() - h->next()->vertex()->point());		

		// the plane defined by the two edges

		vec3 plane_n = normalize(cross(plane_v1, plane_v0));


		// orthogonalize in-plane vectors

		plane_v0 = normalize(plane_v0 - dot(plane_v0, plane_v1) * plane_v1);

		vec3 plane_origin = h->next()->vertex()->point();


		// 2d coordinates in plane

		vec3 local_t = local_plane_coords(h->vertex()->point(), plane_v0, plane_v1, plane_n, plane_origin);
		vec3 local_b = local_plane_coords(h->opposite()->vertex()->point(), plane_v0, plane_v1, plane_n, plane_origin);
		vec3 local_l = local_plane_coords(h->next()->vertex()->point(), plane_v0, plane_v1, plane_n, plane_origin);
		vec3 local_r = local_plane_coords(h->opposite()->next()->vertex()->point(), plane_v0, plane_v1, plane_n, plane_origin);


		// check if edge intersections lies inside triangles pair (in plane)

		vec3 tb = local_t - local_b;
		vec3 lr = local_l - local_r;

		double ntb[2];
		double ctb;

		double nlr[2];
		double clr;

		ntb[0] = - tb[1];
		ntb[1] = tb[0];

		ctb = -(ntb[0] * local_t[0] + ntb[1] * local_t[1]);


		nlr[0] = - lr[1];
		nlr[1] = lr[0];

		clr = -(nlr[0] * local_l[0] + nlr[1] * local_l[1]);


		double det = ntb[0] * nlr[1] - nlr[0] * ntb[1];

		vec3 intersection(- (nlr[1] * ctb - ntb[1] * clr) / det, 
			- (-nlr[0] * ctb + ntb[0] * clr) / det, 
			0.0);


		double l0 = dot(intersection - local_r, lr) / dot(lr, lr);
		double l1 = dot(intersection - local_b, tb) / dot(tb, tb);

		return l0 > 0.0 && l0 < 1.0 && l1 > 0.0 && l1 < 1.0;
	}

	vec3 IrregEditor::local_plane_coords(
		const vec3& p, const vec3& v0, const vec3& v1, const vec3& plane_n, const vec3& plane_origin
		) {
			vec3 q(p + dot(plane_origin - p, plane_n) * plane_n);
			return vec3( dot(q - plane_origin, v0), dot(q - plane_origin, v1), 0.0);
	}

	void IrregEditor::split_vertex(Vertex* v, Halfedge* f1, Halfedge* g1) {
		Vertex *v2 = new_vertex() ;
		Halfedge *f2 = new_edge() ;
		Halfedge *g2 = new_edge() ;
		Halfedge *f =  new_edge() ;
		Halfedge *g =  f->opposite() ;
		std::vector<Halfedge*> fan ;
		//dxy test
		std::cout << "in split_vertex():" << std::endl;
		std::cout << "v: (" << v->point().x << ", " << v->point().y << ")" << std::endl;
		std::cout << "f1: (" << f1->vertex()->point().x << ", " << f1->vertex()->point().y << ") --> (";
		std::cout << f1->opposite()->vertex()->point().x << ", " << f1->opposite()->vertex()->point().y << ")" << std::endl;
		std::cout << "g1: (" << g1->vertex()->point().x << ", " << g1->vertex()->point().y << ") --> (";
		std::cout << g1->opposite()->vertex()->point().x << ", " << g1->opposite()->vertex()->point().y << ")" << std::endl;
		//
		gx_assert(f1->opposite()->vertex()==v && g1->opposite()->vertex()==v) ;

		Halfedge *cir = g1->opposite() ;
		do {
			fan.push_back(cir) ;
			cir = cir->opposite()->prev() ;
		} while (cir != f1->opposite()) ;


		link(f2, f1->opposite()->next(), 1) ;
		link(f1->opposite()->prev(), f2, 1) ;
		link(f2->opposite(), f1->opposite(), 1) ;
		link(f1->opposite(), f, 1) ;
		link(f, f2->opposite(), 1) ;
		set_halfedge_vertex(f2->opposite(), f1->vertex()) ;
		set_halfedge_vertex(f2, v2) ;
		set_halfedge_vertex(f, v2) ;
		gx_assert(f1->opposite()->vertex()==v) ;
		set_facet_on_orbit(f, new_facet()) ;
		make_facet_key(f, f->facet()) ;
		make_vertex_key(f) ;
		make_vertex_key(f2->opposite()) ;

		set_facet_on_orbit(f2, f2->next()->facet()); 
		make_facet_key(f2, f2->facet()) ;

		link(g2, g1->opposite()->next(), 1) ;
		link(g1->opposite()->prev(), g2, 1) ;
		link(g1->opposite(), g, 1) ;
		link(g, g2->opposite(), 1) ;
		link(g2->opposite(), g1->opposite(), 1) ;
		set_halfedge_vertex(g, v) ;
		set_halfedge_vertex(g2, v) ;
		set_halfedge_vertex(g2->opposite(), g1->vertex()) ;
		set_halfedge_vertex(g1->opposite(), v2) ;
		set_facet_on_orbit(g, new_facet()) ;
		make_facet_key(g, g->facet()) ;
		make_vertex_key(g) ;

		set_facet_on_orbit(g2, g2->next()->facet()); 
		make_facet_key(g2, g2->facet()) ;

		for(int i=0; i<fan.size(); ++i) {
			set_halfedge_vertex(fan[i], v2) ;
		}

		///dxy add
		Property<Map::Vertex, bool> is_primary(map_->vertex_property_manager(), "is_primary", PropertyManager::FIND);
		is_primary[v2] = is_primary[v];
		///

		// smooth
		v2->set_point(v->point()) ;
		std::vector<Vertex*> to_smth ;
		to_smth.push_back(v) ;
		to_smth.push_back(v2) ;
		smooth_vertices(to_smth, 5) ;
	}

	void IrregEditor::smooth(int iter) {
		std::vector<vec3> new_points ;

		for(int i=0; i<iter; ++i) {
			int cur = 0 ;

			new_points.assign(map_->nb_vertices(), vec3(0, 0, 0)) ;
			FOR_EACH_VERTEX(Map, map_, v) {
				Halfedge *h = v->halfedge() ;
				if(!v->is_on_border()) {
					do {
						new_points[cur] += h->opposite()->vertex()->point() ;
						h = h->next_around_vertex() ;
					} while(h!=v->halfedge()) ;
					new_points[cur] = 1.0/v->degree()*new_points[cur] ;
				}
				else {
					new_points[cur] = v->point() ;
				}
				cur ++ ;
			}

			cur = 0 ;
			FOR_EACH_VERTEX(Map, map_, v) {
				v->set_point(new_points[cur]) ;
				cur ++ ;
			}
		}

	}

	void IrregEditor::smooth_vertices(std::vector<Vertex*>& to_smth, int iter) {
		std::vector<vec3> new_points ;
		new_points.resize(to_smth.size()) ;
		for(int i=0; i<iter; ++i) {
			new_points.assign(new_points.size(), vec3(0, 0, 0)) ;

			for(unsigned int j=0; j<to_smth.size(); ++j) {
				Halfedge *h = to_smth[j]->halfedge() ;
				do {
					new_points[j] += h->opposite()->vertex()->point() ;
					h = h->next_around_vertex() ;
				} while(h!=to_smth[j]->halfedge()) ;
				new_points[j] = 1.0/to_smth[j]->degree()*new_points[j] ;
			}

			for(unsigned int j=0; j<to_smth.size(); ++j) {
				to_smth[j]->set_point(new_points[j]) ;
			}
		}
	}

	bool IrregEditor::inside_facet(Facet* f, vec2& pt) {
		Halfedge* h = f->halfedge() ;
		vec3 v0 = h->vertex()->point() ;
		vec3 v1 = h->next()->vertex()->point() ;
		vec3 v2 = h->prev()->vertex()->point() ;
		vec3 p(pt.x, pt.y, 0) ;

		vec3 c01 = cross(v0-p, v1-p) ;
		vec3 c12 = cross(v1-p, v2-p) ;
		vec3 c20 = cross(v2-p, v0-p) ;

		if(dot(c01, c12)>0 && dot(c01, c20)>0) {
			return true ;
		}
		return false ;
	}

	Map::Facet* IrregEditor::locate(vec2& pt) {
		FOR_EACH_FACET(Map, map_, it) {
			if(inside_facet(it, pt)) {
				return it ;
			}
		}
		return nil ;
	}

	void IrregEditor::pick_vertex(vec2& pt) {
		Facet* f = locate(pt) ;
		cur_f_ = f ;
		vec3 p3d(pt.x, pt.y, 0) ;
		double min_dist = 1e10 ;
		Vertex* v = nil ;
		if(f!=nil) {
			Halfedge* h = f->halfedge() ;
			do {
				double cur_dist = distance(h->vertex()->point(), p3d) ;
				if(cur_dist < min_dist && cur_dist < 1e-2) {
					min_dist = cur_dist ;
					v = h->vertex() ;
				}
				h = h->next() ;
			} while (h!=f->halfedge()) ;
		}

		if(v!=nil && std::find(picked_verts_.begin(), picked_verts_.end(), v)==picked_verts_.end()) {
			picked_verts_.push_back(v) ;
		}
	}

	static inline double distance_point_segment(vec3& v0, vec3& v1, vec3& p) {
		if(dot(p-v0, v1-v0) > 0) {
			return cross(normalize(v1-v0), p-v0).length() ;
		} else {
			return gx_min(distance(p, v0), distance(p, v1)) ;
		}
	}

	void IrregEditor::pick_edge(vec2& pt) {
		Facet* f = locate(pt) ;
		cur_f_ = f ;
		if(f==nil) return ;
		vec3 p3d(pt.x, pt.y, 0) ;
		double min_dist = 1e10 ;
		Halfedge *e = nil ;

		Halfedge* h = f->halfedge() ;
		do {
			double cur_dist = distance_point_segment(h->vertex()->point(), h->opposite()->vertex()->point(), p3d) ;
			if(cur_dist < min_dist) {
				min_dist = cur_dist ;
				e = h ;
			}
			h = h->next() ;
		} while (h!=f->halfedge()) ;

		if(e!=nil) {
			bool flag1 = std::find(picked_edges_.begin(), picked_edges_.end(), e)==picked_edges_.end() ;
			bool flag2 = std::find(picked_edges_.begin(), picked_edges_.end(), e->opposite())==picked_edges_.end() ;
			if(flag1 && flag2) {
				picked_edges_.push_back(e) ;
			}
		}
	}

	static inline int vindex_grid(int i, int j, int n, int margin) {
		return (n+2*margin)*(i+margin) + j+margin ;
	}

	void IrregEditor::create_grid_map(double radius) {
		int    n = ceil(1.0/radius) ;
		double len = 1.0/n ;
		int    margin = 3 ;

		map_->clear() ;
		MapBuilder builder(map_) ;
		builder.begin_surface() ;

		for(int i=-margin; i<n+margin; ++i) {
			for(int j=-margin; j<n+margin; ++j) {
				builder.add_vertex(vec3((i+0.5)*len, (j+0.5)*len, 0)) ;
			}
		}

		for(int i=-margin; i<n+margin-1; ++i) {
			for(int j=-margin; j<n+margin-1; ++j) {
				builder.begin_facet() ;
				builder.add_vertex_to_facet(vindex_grid(i  , j  , n, margin)) ;
				builder.add_vertex_to_facet(vindex_grid(i+1, j  , n, margin)) ;
				builder.add_vertex_to_facet(vindex_grid(i+1, j+1, n, margin)) ;
				builder.end_facet() ;

				builder.begin_facet() ;
				builder.add_vertex_to_facet(vindex_grid(i  , j  , n, margin)) ;
				builder.add_vertex_to_facet(vindex_grid(i+1, j+1, n, margin)) ;
				builder.add_vertex_to_facet(vindex_grid(i  , j+1, n, margin)) ;
				builder.end_facet() ;
			}
		}

		builder.end_surface() ;	
	}

	static inline int vindex_regular(int i, int j, int n, int m, int margin) {
		return (n+2*margin)*(i+margin) + j+margin ;
	}

	void IrregEditor::create_regular_map(double radius) {
		int n = 1.0 / radius ;
		int m = int(2.0/sqrt(3.0)*n+0.5) ;
		double newr = sqrt(2.0/(sqrt(3.0)*m*n)) ;
		double width = n*newr ;
		double height = m*sqrt(3.0)/2*newr ;
		int margin = 3 ;

		map_->clear() ;
		MapBuilder builder(map_) ;
		builder.begin_surface() ;

		for(int i=-margin; i<m+margin; ++i) {
			for(int j=-margin; j<n+margin; ++j) {
				if(i%2==0) {
					builder.add_vertex(vec3(j*newr, i*sqrt(3.0)/2*newr, 0)) ;
				} 
				else {
					builder.add_vertex(vec3((j+0.5)*newr, i*sqrt(3.0)/2*newr, 0)) ;
				}				
			}
		}

		for(int i=-margin; i<m+margin-1; ++i) {
			for(int j=-margin; j<n+margin-1; ++j) {
				builder.begin_facet() ;
				if(i%2==0) {
					builder.add_vertex_to_facet(vindex_regular(i  , j  , n, m, margin)) ;
					builder.add_vertex_to_facet(vindex_regular(i  , j+1, n, m, margin)) ;
					builder.add_vertex_to_facet(vindex_regular(i+1, j  , n, m, margin)) ;
				}else {
					builder.add_vertex_to_facet(vindex_regular(i  , j  , n, m, margin)) ;
					builder.add_vertex_to_facet(vindex_regular(i  , j+1, n, m, margin)) ;
					builder.add_vertex_to_facet(vindex_regular(i+1, j+1, n, m, margin)) ;
				}
				builder.end_facet() ;

				builder.begin_facet() ;
				if(i%2==0) {
					builder.add_vertex_to_facet(vindex_regular(i  , j+1, n, m, margin)) ;
					builder.add_vertex_to_facet(vindex_regular(i+1, j+1, n, m, margin)) ;
					builder.add_vertex_to_facet(vindex_regular(i+1, j  , n, m, margin)) ;
				}else {
					builder.add_vertex_to_facet(vindex_regular(i  , j  , n, m, margin)) ;
					builder.add_vertex_to_facet(vindex_regular(i+1, j+1, n, m, margin)) ;
					builder.add_vertex_to_facet(vindex_regular(i+1, j  , n, m, margin)) ;
				}
				builder.end_facet() ;
			}
		}

		builder.end_surface() ;	
	}

	void IrregEditor::perturb(double esp) {
		FOR_EACH_VERTEX(Map, map_, it) {
			vec3 dir = normalize(vec3(Numeric::random_float64(), Numeric::random_float64(), 0)) ;
			it->set_point(it->point() + esp*dir) ;
		}
	}

	///dxy add
	Point IrregEditor::linear_mix(const Point &s, const Point &t, double q) {
		//q=0, return s; 
		//q=1, return t;
		return to_cgal(q * to_geex(t) + (1-q) * to_geex(s));
	}

	int IrregEditor::delaunay_status(const std::vector<Point> &pts) {
		gx_assert(pts.size()==4);
		// convex test
		for(int i=0; i<4; ++i) {
			if (!CGAL::left_turn(pts[i], pts[(i+1)%4], pts[(i+2)%4])) {
				std::cout <<"pts:" << std::endl;
				for(int j=0; j<4; ++j) {
					std::cout << "(" << pts[j].x() << "," << pts[j].y() << ") ";
				}
				std::cout << std::endl;
				std::cout << "non-left turn:" << i << std::endl;
				
				return -1;  //non-convex
			}
		}
		// delaunay test
		if ((CGAL::side_of_oriented_circle(pts[0], pts[1], pts[2], pts[3]) == CGAL::ON_NEGATIVE_SIDE)
			&& (CGAL::side_of_oriented_circle(pts[2], pts[3], pts[0], pts[1]) == CGAL::ON_NEGATIVE_SIDE)) {
				return 1;  //convex and delaunay
		} else return 0;   //convex and non-delaunay
	}

	int IrregEditor::delaunay_status(const std::vector<Point> &source, const std::vector<Point> &target, double q) {
		gx_assert(source.size()==4);
		gx_assert(target.size()==4);
		std::vector<Point> pts;
		//
// 		std::cout <<"src:" << std::endl;
// 		for(int j=0; j<4; ++j) {
// 			std::cout << "(" << source[j].x() << "," << source[j].y() << ") ";
// 		}
// 		std::cout << std::endl;
// 		//
// 		std::cout <<"dst:" << std::endl;
// 		for(int j=0; j<4; ++j) {
// 			std::cout << "(" << target[j].x() << "," << target[j].y() << ") ";
// 		}
// 		std::cout << std::endl;
		//
		//std::cout << "q: " << q << std::endl;
		//
		for(int i=0; i<4; ++i) {
			pts.push_back(linear_mix(source[i], target[i], q));
		}
		return delaunay_status(pts);
	}

	bool IrregEditor::is_delaunay_edge(Map::Halfedge_iterator eit) {
		if (eit->is_border_edge()) return true;
		else {
			std::vector<Point> pts;
			pts.push_back(to_cgal(vec2(eit->vertex()->point().x, eit->vertex()->point().y)));
			pts.push_back(to_cgal(vec2(eit->next()->vertex()->point().x, eit->next()->vertex()->point().y)));
			pts.push_back(to_cgal(vec2(eit->prev()->vertex()->point().x, eit->prev()->vertex()->point().y)));
			pts.push_back(to_cgal(vec2(eit->opposite()->next()->vertex()->point().x, eit->opposite()->next()->vertex()->point().y)));

			return ((CGAL::side_of_oriented_circle(pts[0], pts[1], pts[2], pts[3]) == CGAL::ON_NEGATIVE_SIDE)
				&& (CGAL::side_of_oriented_circle(pts[2], pts[3], pts[0], pts[1]) == CGAL::ON_NEGATIVE_SIDE));
			//return (delaunay_status(pts) == 1);
		}
	}

	std::vector<Map::Halfedge_iterator> IrregEditor::collect_non_delaunay_edges() {
		return std::vector<Map::Halfedge_iterator>(); //TODO: implement
	}

	void IrregEditor::topo_delaunay_optimize() {
		//assert: periodic triangulation, has_vertex_property(is_primary, primary_id)
		std::vector<Map::Halfedge_iterator> vNDE = collect_non_delaunay_edges();
		//std::set<Map::Halfedge_iterator> NDE(vNDE.begin(), vNDE.end());  //the less object must be defined
		while (vNDE.size()>0)
		{
			//TODO: delaunay optimize

		}
	}


	//Naive Version
	void IrregEditor::topo_delaunay_optimize(unsigned int nb_max_iter) {
		std::cout << "topo delaunay optimize begin..." << std::endl;
		Property<Map::Vertex, bool> is_primary(map_->vertex_property_manager(), "is_primary", PropertyManager::FIND);
		Property<Map::Vertex, int> primary_id(map_->vertex_property_manager(), "primary_id", PropertyManager::FIND);
		std::cout << "vertex property found." << std::endl;
		//build a map from primary_id to the set of vertices with that id
		std::map<int, std::set<Map::Vertex*>> id2vertices;
		FOR_EACH_VERTEX(Map, map_, v) {
			int pid = primary_id[v];
			if (id2vertices.find(pid) == id2vertices.end()) {
				id2vertices.insert(std::make_pair(pid, std::set<Map::Vertex*>()));
			} 
			id2vertices[pid].insert(v);
		}
		std::cout << "id2vertices built" << std::endl;
		//delaunay optimize
		unsigned int nb_iter = 0;
		int nb_non_delaunay_edges = 0;
		do 
		{
			nb_non_delaunay_edges = 0;
			FOR_EACH_EDGE(Map, map_, eit) {
				if (eit->is_border_edge()) continue;
				Map::Vertex* v0 = eit->vertex();
				Map::Vertex* v1 = eit->opposite()->vertex();
				if((is_primary[v0] && is_primary[v1])
					|| (is_primary[v0] && primary_id[v0]<primary_id[v1])
					|| (is_primary[v1] && primary_id[v1]<primary_id[v0]))
				{ //eit is primary edge
					if (!is_delaunay_edge(eit)) {
						std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
						std::cout << "a primary non-delaunay edge found" << std::endl;
						++nb_non_delaunay_edges;
						//
						std::vector<Map::Vertex*> vs;
						vs.push_back(eit->vertex());
						vs.push_back(eit->next()->vertex());
						vs.push_back(eit->prev()->vertex());
						vs.push_back(eit->opposite()->next()->vertex());
						//
						std::vector<Point> pts;
						for(int i=0; i<4; ++i) {
							pts.push_back(to_cgal(vec2(vs[i]->point().x, vs[i]->point().y)));
						}
						//assert: pts[0][1][2][3] is a convex quad
						gx_assert(delaunay_status(pts) != -1);
						//compute target points
						std::vector<Circle_2> circles;
						std::vector<Point> circle_centers;
						std::vector<double> circle_radius2;
						for(int i=0; i<4; ++i) {
							circles.push_back(Circle_2(pts[(i+1)%4], pts[(i+2)%4], pts[(i+3)%4]));
							circle_centers.push_back(circles.back().center());
							circle_radius2.push_back(circles.back().squared_radius());
						}
						//
						std::vector<Point> target_pts = pts;
						//type1 target pts
// 						if (circles[0].bounded_side(pts[0]) == CGAL::ON_UNBOUNDED_SIDE) {
// 							//manually compute the intersection between circles[0] and ray(circles[0].center(), pts[0])
// 							double q = sqrt(circle_radius2[0] / distance2(to_geex(circle_centers[0]), to_geex(pts[0])));
// 							target_pts[0] = linear_mix(circle_centers[0], pts[0], q);
// 						}
// 						if (circles[1].bounded_side(pts[1]) == CGAL::ON_BOUNDED_SIDE) {
// 							double q = sqrt(circle_radius2[1] / distance2(to_geex(circle_centers[1]), to_geex(pts[1])));
// 							target_pts[1] = linear_mix(circle_centers[1], pts[1], q);
// 						}
// 						if (circles[2].bounded_side(pts[2]) == CGAL::ON_UNBOUNDED_SIDE) {
// 							double q = sqrt(circle_radius2[2] / distance2(to_geex(circle_centers[2]), to_geex(pts[2])));
// 							target_pts[2] = linear_mix(circle_centers[2], pts[2], q);
// 						}
// 						if (circles[3].bounded_side(pts[3]) == CGAL::ON_BOUNDED_SIDE) {
// 							double q = sqrt(circle_radius2[3] / distance2(to_geex(circle_centers[3]), to_geex(pts[3])));
// 							target_pts[3] = linear_mix(circle_centers[3], pts[3], q);
// 						}
						//type2 target pts
						for(int i=0; i<4; ++i) {
							vec2 Pi = to_geex(pts[i]);
							vec2 Pk = to_geex(pts[(i+2)%4]);
							vec2 Oi = to_geex(circle_centers[i]);
							vec2 PkOi = Oi - Pk;
							vec2 PkPi = Pi - Pk;
							vec2 target_i = Pk + 2 * PkPi * dot(PkOi, PkPi) / dot(PkPi, PkPi);
							target_pts[i] = to_cgal(target_i);
						}
						//print src and dst points
						std::cout <<"src = {" ;
						for(int i=0; i<3; ++i) {
							std::cout << "(" << pts[i].x() << "," << pts[i].y() << ") ,";
						}
						std::cout << "(" << pts[3].x() << "," << pts[3].y() << ")}" << std::endl;
						std::cout <<"dst = {" ;
						for(int i=0; i<3; ++i) {
							std::cout << "(" << target_pts[i].x() << "," << target_pts[i].y() << ") ,";
						}
						std::cout << "(" << target_pts[3].x() << "," << target_pts[3].y() << ")}" << std::endl;
						// search for delaunay 
						std::cout << "search for delaunay..." << std::endl;
						double lo = 0;
						double hi = 1;
						double mi = hi;
						while (true) {
							//Sleep(200);
							std::cout << "mi: " << mi << std::endl;
							bool found = false;
							switch (delaunay_status(pts, target_pts, mi))
							{
							case -1: //non-convex
								std::cout << "non-convex" << std::endl;
								hi = mi;
								mi = 0.5 * (hi + lo);
								break;
							case 0: //non-delaunay
								std::cout << "non-delaunay" << std::endl;
								lo = mi;
								mi = 0.5 * (hi + lo);
								break;
							case 1: //delaunay
								found = true;
								break;
							default:
								std::cerr << "topo_delaunay_optimize::Error: invalid delaunay status" << std::endl;
								return;
							}
							if (found) break;
						}
						double q_upper = mi;
						std::cout << "q_upper: " << q_upper << std::endl; 
						// search for minimal delaunay
						std::cout << "search for minimal delaunay..." << std::endl;
						hi = mi;
						mi = 0.5 * (hi + lo);
						while(hi-lo > 0.01) {
							//Sleep(200);
							std::cout << "(" << lo << ", " << hi << "]" << std::endl;
							switch (delaunay_status(pts, target_pts, mi))
							{
							case 0: //non-delaunay
								lo = mi;
								mi = 0.5 * (hi + lo);
								break;
							case 1: //delaunay
								hi = mi;
								mi = 0.5 * (hi + lo);
								break;
							default:
								std::cerr << "topo_delaunay_optimize::Error: unexpected delaunay status" << std::endl;
								return;
							}
						}
						double q_lower = hi;
						std::cout << "q_lower: " << q_lower << std::endl;
						std::cout << "q -> [" << q_lower << ", " << q_upper << "]" << std::endl;
						gx_assert(q_upper >= q_lower); //dxy test
						//modify q
						double q_star = 0.9 * q_lower + 0.1 * q_upper;
						std::cout << "q*: " << q_star << std::endl;
						//print the result pts
						std::cout <<"res = {" ;
						for(int i=0; i<3; ++i) {
							std::cout << "(" << (1-q_star) * pts[i].x() + q_star * target_pts[i].x() << "," 
								<< (1-q_star) * pts[i].y() + q_star * target_pts[i].y() << ") ,";
						}
						std::cout << "(" << (1-q_star) * pts[3].x() + q_star * target_pts[3].x() << "," 
							<< (1-q_star) * pts[3].y() + q_star * target_pts[3].y() << ")}" << std::endl;
						std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
						//(periodically) apply vertex movement according to q
						std::vector<vec3> displacement;
						for(int i=0; i<4; ++i) {
							vec2 dis = q_star * (to_geex(target_pts[i]) - to_geex(pts[i]));
							displacement.push_back(vec3(dis.x, dis.y, 0));
						}
						int num_displaced_points = 0;
						for(int i=0; i<4; ++i) {
							const std::set<Map::Vertex*> &period_verts = id2vertices[primary_id[vs[i]]];
							for(auto it=period_verts.begin(); it!=period_verts.end(); ++it) {
								(*it)->set_point((*it)->point() + displacement[i]);
								++num_displaced_points;
							}
						}
						std::cout << num_displaced_points << " points displaced." << std::endl;
					}
				}
			}
			std::cout << "Iter " << nb_iter << ": " << nb_non_delaunay_edges << " NDE found." << std::endl;
			++nb_iter;
		} while (nb_non_delaunay_edges>0 && nb_iter<nb_max_iter);
	}

	
	void IrregEditor::do_editing_from_delaunay(const std::vector<int> &vertices, const std::vector<std::pair<int,int>> &edges, int nb_max_iter /*=50*/) {
		delaunay_to_map();
		//
		Property<Map::Vertex, bool> is_primary(map_->vertex_property_manager(), "is_primary", PropertyManager::FIND);
		Property<Map::Vertex, int> primary_id(map_->vertex_property_manager(), "primary_id", PropertyManager::FIND);
		//build a map from primary_id to the set of vertices with that id
		std::map<int, std::set<Map::Vertex*>> id2vertices;
		FOR_EACH_VERTEX(Map, map_, v) {
			int pid = primary_id[v];
			if (id2vertices.find(pid) == id2vertices.end()) {
				id2vertices.insert(std::make_pair(pid, std::set<Map::Vertex*>()));
			} 
			id2vertices[pid].insert(v);
		}
		//
		clear_picked();
		//pick vertices
		for (auto it=vertices.begin(); it!=vertices.end(); ++it) {
			int vid = *it;
			std::set<Map::Vertex*> vs = id2vertices[vid];
			for(auto vit=vs.begin(); vit!=vs.end(); ++vit) {
				picked_verts_.push_back(*vit);
			}
		}
		//pick edges
		for (auto eit=edges.begin(); eit!=edges.end(); ++eit) {
			int id0 = eit->first;
			int id1 = eit->second;
			std::set<Map::Vertex*> v0s = id2vertices[id0];
			for(auto vit=v0s.begin(); vit!=v0s.end(); ++vit) {
				Map::Halfedge *h = (*vit)->halfedge() ;
				do {
					if(primary_id[h->opposite()->vertex()] == id1) {
						picked_edges_.push_back(h);
						break;
					}
					h = h->next_around_vertex() ;
				} while(h!=(*vit)->halfedge()) ;
			}
		}
		//do editing
		std::cout << "before editing, " << map_->nb_vertices() << " map vertices." << std::endl;
		do_editing();
		//topo delaunay optimize
		std::cout << "after editing, " << map_->nb_vertices() << " map vertices." << std::endl;
		topo_delaunay_optimize(nb_max_iter);
		std::cout << "after TDO, " << map_->nb_vertices() << " map vertices." << std::endl;
		//
		map_to_delaunay();
	}

	void IrregEditor::flip_edges(const std::vector<std::pair<int,int>> &edges, int nb_max_iter /*=50*/) {
		delaunay_to_map();
		//
		Property<Map::Vertex, bool> is_primary(map_->vertex_property_manager(), "is_primary", PropertyManager::FIND);
		Property<Map::Vertex, int> primary_id(map_->vertex_property_manager(), "primary_id", PropertyManager::FIND);
		//build a map from primary_id to the set of vertices with that id
		std::map<int, std::set<Map::Vertex*>> id2vertices;
		FOR_EACH_VERTEX(Map, map_, v) {
			int pid = primary_id[v];
			if (id2vertices.find(pid) == id2vertices.end()) {
				id2vertices.insert(std::make_pair(pid, std::set<Map::Vertex*>()));
			} 
			id2vertices[pid].insert(v);
		}
		//
		std::vector<Map::Halfedge*> map_edges;
		for (auto eit=edges.begin(); eit!=edges.end(); ++eit) {
			int id0 = eit->first;
			int id1 = eit->second;
			std::set<Map::Vertex*> v0s = id2vertices[id0];
			for(auto vit=v0s.begin(); vit!=v0s.end(); ++vit) {
				Map::Halfedge *h = (*vit)->halfedge() ;
				do {
					if(primary_id[h->opposite()->vertex()] == id1) {
						map_edges.push_back(h);
						break;
					}
					h = h->next_around_vertex() ;
				} while(h!=(*vit)->halfedge()) ;
			}
		}
		std::cout << map_edges.size() << " map edges to flip" << std::endl;
		//
		for(auto eit=map_edges.begin(); eit!=map_edges.end(); ++eit) {
			edge_flip(*eit);
		}
		std::cout << "edges flipped" << std::endl;
		topo_delaunay_optimize(nb_max_iter);
		//
		map_to_delaunay();
	}

	void IrregEditor::test() {
		delaunay_to_map();
		//
		Property<Map::Vertex, bool> is_primary(map_->vertex_property_manager(), "is_primary", PropertyManager::FIND);
		Property<Map::Vertex, int> primary_id(map_->vertex_property_manager(), "primary_id", PropertyManager::FIND);
		//build a map from primary_id to the set of vertices with that id
		std::map<int, std::set<Map::Vertex*>> id2vertices;
		FOR_EACH_VERTEX(Map, map_, v) {
			int pid = primary_id[v];
			if (id2vertices.find(pid) == id2vertices.end()) {
				id2vertices.insert(std::make_pair(pid, std::set<Map::Vertex*>()));
			} 
			id2vertices[pid].insert(v);
		}
		//
		Map::Halfedge_iterator e;
		FOR_EACH_EDGE(Map, map_, eit) {
			if (eit->is_border_edge()) continue;
			Map::Vertex* v0 = eit->vertex();
			Map::Vertex* v1 = eit->opposite()->vertex();
			if((is_primary[v0] && is_primary[v1])
				|| (is_primary[v0] && primary_id[v0]<primary_id[v1])
				|| (is_primary[v1] && primary_id[v1]<primary_id[v0]))
			{
				e = eit;
				break;
			}
		}
		//collect mirrors of e
		int pid0 = primary_id[e->vertex()];
		int pid1 = primary_id[e->opposite()->vertex()];
		std::vector<Map::Halfedge_iterator> all_e;
		FOR_EACH_HALFEDGE(Map, map_, eit) {
			if (eit->is_border_edge()) continue;
			Map::Vertex* v0 = eit->vertex();
			Map::Vertex* v1 = eit->opposite()->vertex();
			if ((primary_id[v0]==pid0) && (primary_id[v1]==pid1)) {
				all_e.push_back(eit);
			}
		}
		std::cout << "found " << all_e.size() << " copies of e." << std::endl;
		//
		for(auto it=all_e.begin(); it!=all_e.end(); ++it) {
			edge_flip(*it);
		}
		std::cout << "edges flipped" << std::endl;
		topo_delaunay_optimize(50);
		//
		for(auto it=all_e.begin(); it!=all_e.end(); ++it) {
			edge_flip(*it);
		}
		std::cout << "flipped back" << std::endl;
		topo_delaunay_optimize(10);
	}
	///dxy add end

} // end of namespace Geex		