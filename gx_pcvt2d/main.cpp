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

#include <Geex/graphics/geexapp.h>
#include <Geex/basics/file_system.h>
#include <Geex/basics/stopwatch.h>
#include <glut_viewer/tweak_bar.h>
#include <fstream>
#include "cvt.h"
//#include "psa-src\psa.h"
//#include "irreg_editor.h"

void Lloyd() ;
void TW_CALL CB_action(void *clientData) ;
void TW_CALL CB_convert_mesh(void *clientData) ;
void TW_CALL CB_create_grid(void *clientData) ;

namespace Geex {

	class CVTApp : public GeexApp {
    public:
        CVTApp(int argc, char** argv) : GeexApp(argc, argv) { 
            hdr_ = false ;
            boundary_filename_ = get_file_arg("line") ;
            if(boundary_filename_.length() > 0) {
                if(!Geex::FileSystem::is_file(boundary_filename_)) {
                    boundary_filename_ = Geex::FileSystem::get_project_root() + 
                        "/gx_pcvt2d/" + boundary_filename_ ;
                }
            }
            nb_points_ = 100 ;
            get_arg("nb_pts", nb_points_) ;
            nb_iter_ = 50 ;
            get_arg("nb_iter", nb_iter_) ;            
            non_convex_ = GL_FALSE ;
            get_arg("non_convex", non_convex_) ;
            edit_ = GL_FALSE ;
			///dxy add
			select_ = GL_FALSE;
			///
			insert_boundary_ = GL_FALSE ;
			get_arg("insert_boundary", insert_boundary_) ;
			min_max_ = 30 ;
			nb_runs_ = 1 ;
        }

        CVT* cvt() { return static_cast<CVT*>(scene()) ; }

        GLboolean& edit() { return edit_ ; }
		///dxy
		GLboolean& select() { return select_ ; }

		virtual void init_scene() {
            scene_ = new CVT ;
            std::cerr << "Non convex = " 
                      << (non_convex_ ? "true" : "false") 
                      << "(use +non_convex to set)" << std::endl ;
            cvt()->set_non_convex_mode(non_convex_) ;
            if(boundary_filename_.length() > 0) {
                cvt()->load_boundary(boundary_filename_) ;
            }
			cvt()->insert_boundary() = insert_boundary_ ;
            cvt()->insert_random_vertices(nb_points_) ;
			insert_copies() ;
        }

		void insert_grid() {
			cvt()->clear() ;
			cvt()->insert_grid(500) ;
			insert_copies() ;
		}

        void Lloyd() {
            cvt()->lloyd(nb_iter_) ;
        }

        void Lloyd(int nb_iter) {
            cvt()->lloyd(nb_iter) ;
        }

		///dxy add
		void Update_Energy_Grad() {
			std::cout << "Energy = " << cvt()->lloyd_energy() << std::endl ;
		}

		void Update_Total_Area() {
			cvt()->delaunay()->update_area();
		}

		GLfloat& nb_iter() { return nb_iter_ ; }

		void NewtonLloyd() {
			cvt()->lloyd_energy(); //trick: otherwise, update_area() will not work
			cvt()->delaunay()->update_area();
			if (cvt()->use_topo_optimization())
			{
				cvt()->newton_lloyd(1);
				int it = 0;
				while (it < nb_iter_)
				{
					if (cvt()->period())
					{
						if (cvt()->use_balanced_stretch())
						{
							cvt()->stretch_topo_optimize_period_balanced();
						}
						else {
							cvt()->stretch_topo_optimize_period();
						}
					}
					else
					{
						cvt()->stretch_topo_optimize();
					}
					cvt()->newton_lloyd(5);
					//test info
					cvt()->print_topo_info(it);
					cvt()->calc_separatrix_graph();
					it += 1;
				}
			}
			else {
				cvt()->newton_lloyd(nb_iter_) ;
			}
		}

		void pair57_move() {
			std::vector<Pair57> P;
			cvt()->collect_selected_Pair57(P);
			if (!P.empty())
			{
				std::vector<PairMoveCode> M;
				cvt()->collect_Pair57_operation(P, M);
				for (int i=0; i<cvt()->num_Pair57_move(); ++i)
				{
					std::cout << "5-7 Pair Move : " << i << std::endl;
					cvt()->move_Pair57(P, M);
					cvt()->newton_lloyd(cvt()->num_Pair57_smoothing()); //smoothing
					int num_valid = 0;
					for (int j=0; j<P.size(); ++j)
					{
						if (P[j].verify())
						{
							num_valid++;
						}
					}
					std::cout << num_valid << " valid pairs left"<< std::endl;
					if (num_valid == 0) break;
				}
			}
		}
		///

//         void NewtonLloyd() {
//             cvt()->newton_lloyd(nb_iter_) ;
//         }

		void Lloyd_fpo() {
			cvt()->lloyd_fpo(nb_iter_) ;
		}

        void reset() {
            cvt()->clear() ;
            cvt()->insert_random_vertices(nb_points_) ;
			insert_copies() ; 
        }

		void insert_copies() {
			cvt()->insert_copies(cvt()->pvd_mode()==FULL_COPY) ;
		}

		void clear_copies() {
			cvt()->clear_copies(true) ;
		}

		void save() {
			std::string filename = Geex::FileSystem::get_project_root() + 
				"/gx_pcvt2d/" + "out." ;
			//dxy add
			std::cout << "please input an identifier:";
			std::string file_identifier;
			std::cin >> file_identifier;
			filename = filename + file_identifier + ".pts";
			std::cout << "saving to: " << filename << std::endl;
			//
			cvt()->save(filename) ;
		}
		void load() {
			std::string filename = Geex::FileSystem::get_project_root() + 
				"/gx_pcvt2d/" + "out." ;
			//dxy add
			std::cout << "please input the identifier:";
			std::string file_identifier;
			std::cin >> file_identifier;
			filename = filename + file_identifier + ".pts";
			//
			cvt()->load(filename) ;
		}

		void save(const char* name) {
			std::string filename = Geex::FileSystem::get_project_root() + 
				"/gx_pcvt2d/" + name + ".pts" ;
			cvt()->save(filename) ;
		}
		void load(const char* name) {
			std::string filename = Geex::FileSystem::get_project_root() + 
				"/gx_pcvt2d/" + name + ".pts" ;
			cvt()->load(filename) ;
		}

		void generate_poisson_disk()  {
			std::string filename = Geex::FileSystem::get_project_root() + 
				"/gx_pcvt2d/pointsets/random/" ;
			if(nb_runs_ > 1) 
				cvt()->generate_pointset(cvt()->sample_radius(), nb_runs_, filename) ;
			else {
				cvt()->generate_poisson_disk() ;
			}
		}

		void analysis_pointsets() {
			std::string path = Geex::FileSystem::get_project_root() + 
				"/gx_pcvt2d/pointsets/PVoronoi/r0.005" ;
			cvt()->analyze_pointsets(path) ;
		}

		void smooth() {
			cvt()->smooth(nb_iter_, true) ;
		}
		

		///dxy add: edge topo edit
		void edge_split() {
			if (cvt()->show_long_edge())
			{
				cvt()->split_isolate_long_edge_period();
			}
			else {
				cvt()->split_selected_edge_period();
			}
		}

		void edge_collapse() {
			if (cvt()->show_short_edge())
			{
				cvt()->collapse_isolate_short_edge_period();
			}
			else {
				cvt()->collapse_selected_edge_period();
			}
		}

		void edge_flip() {
			cvt()->flip_selected_edge();
		}
		///dxy add end

        virtual void init_gui() {
            GeexApp::init_gui() ;

            // New-style GUI =====================================================================================

            TwBar* graphics_bar = TwNewBar("Graphics") ;
			TwDefine(" Graphics position='16 10' size='200 480'") ; 
			TwAddVarRW(graphics_bar, "Vertices", TW_TYPE_BOOL8, &cvt()->show_vertices(), "");
            TwAddVarRW(graphics_bar, "Vertex Size", TW_TYPE_FLOAT, &cvt()->vertices_size(), "min=0 max=1 step=0.01") ;
            TwAddVarRW(graphics_bar, "Vertex Color", TW_TYPE_BOOL8, &cvt()->vertices_color(), "") ;
            TwAddVarRW(graphics_bar, "Centers", TW_TYPE_FLOAT, &cvt()->centers_size(), "min=0 max=1 step=0.01") ;
            TwAddVarRW(graphics_bar, "Primal", TW_TYPE_BOOL8, &cvt()->show_primal_mesh(), "") ;
            TwAddVarRW(graphics_bar, "Dual", TW_TYPE_BOOL8, &cvt()->show_dual_mesh(), "") ;
			TwAddVarRW(graphics_bar, "InnerDual", TW_TYPE_BOOL8, &cvt()->show_inner_voronoi(), "") ;
			TwAddVarRW(graphics_bar, "Colorize", TW_TYPE_BOOL8, &cvt()->colorize(), "") ;
			TwAddVarRW(graphics_bar, "Non-hex", TW_TYPE_BOOL8, &cvt()->show_cells(), "") ;
			///dxy add
			TwAddVarRW(graphics_bar, "Mesh", TW_TYPE_BOOL8, &cvt()->show_mesh(), "");
			TwAddVarRW(graphics_bar, "Energy", TW_TYPE_BOOL8, &cvt()->show_energy(), "") ;
			TwAddVarRW(graphics_bar, "Area", TW_TYPE_BOOL8, &cvt()->show_relative_area(), "");  
			//TwAddVarRW(graphics_bar, "Tri Area", TW_TYPE_BOOL8, &cvt()->show_triangle_area(), ""); 
			///
			TwAddVarRW(graphics_bar, "Snap border", TW_TYPE_BOOL8, &cvt()->snap_boundary(), "") ;
			TwAddVarRW(graphics_bar, "Bndcells", TW_TYPE_BOOL8, &cvt()->show_boundary_cells(), "") ;
			TwAddVarRW(graphics_bar, "Euclidean", TW_TYPE_BOOL8, &cvt()->show_pvd_euclidean(), "") ;
			TwAddVarRW(graphics_bar, "Copies", TW_TYPE_BOOL8, &cvt()->show_copies(), "") ;
			TwEnumVal pvd_def[] = { {NO_COPY, "No Copy"}, {FULL_COPY, "Full Copy"}, {MIN_COPY, "Min Copy"}} ;
			TwType tw_pvd = TwDefineEnum("PVDMode", pvd_def, 3) ;
			TwAddVarRW(graphics_bar, "PVD Mode", tw_pvd, &cvt()->pvd_mode(), "") ;
			TwAddVarRW(graphics_bar, "show regularity", TW_TYPE_BOOL8, &cvt()->show_regularity(), "") ;            
            TwAddVarRW(graphics_bar, "Energy", TW_TYPE_BOOL8, &cvt()->show_energy(), "") ;
            TwAddVarRW(graphics_bar, "Bkgnd. field", TW_TYPE_BOOL8, &cvt()->show_field(), "") ;
            TwAddVarRW(graphics_bar, "Quads", TW_TYPE_FLOAT, &cvt()->quad_ratio(), "min=0.0 max=1.5 step=0.01") ;
            TwAddVarRW(graphics_bar, "show Lp", TW_TYPE_BOOL8, &cvt()->show_lp(), "") ;            
            TwAddVarRW(graphics_bar, "Lp view", TW_TYPE_INT32, &cvt()->lp_shader(), "min=0 max=10 step=1") ;
            TwAddVarRW(graphics_bar, "Lp scale", TW_TYPE_FLOAT, &cvt()->new_uniform("scale", 1.0), "min=0.001 max=10 step=0.1") ;
			TwAddVarRW(graphics_bar, "perturb", TW_TYPE_FLOAT, &cvt()->perturb(), "min=0 max=1 step=0.001") ;

			///dxy comment
// 			TwAddSeparator(graphics_bar, "Poisson sampling", "");
// 			TwAddVarRW(graphics_bar, "number of runs", TW_TYPE_INT32, &nb_runs_, "min=1 max=10000 step=1") ;
// 			TwEnumVal pds_def[] = { {PD_CCVT, "ccvt"}, {PD_VORONOI, "Voronoi"}, {PD_BOUNDARY, "Pure"}, {PD_DARTTHROW, "Dart throw"}, {PD_PENROSE, "Penrose"} } ;
//             TwType tw_pds = TwDefineEnum("Sampling mode", pds_def, 5) ;
//             TwAddVarRW(graphics_bar, "Sampling mode", tw_pds, &cvt()->sampling_mode(), "") ;
// 			TwAddVarRW(graphics_bar, "sample radius", TW_TYPE_DOUBLE, &cvt()->sample_radius(), "min=0 max=1 step=0.000001") ;
// 			TwAddVarRW(graphics_bar, "show radius", TW_TYPE_BOOL8, &cvt()->show_disk(), "") ;
// 			TwAddVarRW(graphics_bar, "show min/max", TW_TYPE_BOOL8, &cvt()->show_min_max(), "") ;
			//TwAddVarRW(graphics_bar, "edge histogram", TW_TYPE_BOOL8, &cvt()->show_edge_hist(), "") ; ///dxy comment

			///dxy add
			TwAddSeparator(graphics_bar, "Grad", "");
			TwAddVarRW(graphics_bar, "Lloyd Grad Mag", TW_TYPE_INT32, &cvt()->lloyd_grad_magnification(), "min=1 max=100 step=1");
			TwAddVarRW(graphics_bar, "Direct Grad Mag", TW_TYPE_INT32, &cvt()->direct_grad_magnification(), "min=1 max=100 step=1");
			TwAddVarRW(graphics_bar, "Total Grad Mag", TW_TYPE_INT32, &cvt()->total_grad_magnification(), "min=1 max=100 step=1");
			TwAddVarRW(graphics_bar, "Lloyd Gradient", TW_TYPE_BOOL8, &cvt()->show_lloyd_grad(), "") ;
			TwAddVarRW(graphics_bar, "Direction Gradient", TW_TYPE_BOOL8, &cvt()->show_direction_grad(), "") ;
			TwAddVarRW(graphics_bar, "Total Gradient", TW_TYPE_BOOL8, &cvt()->show_total_grad(), "");
			TwAddVarRW(graphics_bar, "Regular Only", TW_TYPE_BOOL8, &cvt()->show_regular_grad_only(), "");
			TwAddSeparator(graphics_bar, "Edge", "");
			TwAddVarRW(graphics_bar, "Separatrix", TW_TYPE_BOOL8, &cvt()->show_separatrix(), "");
			TwAddVarRW(graphics_bar, "Edge Stretch", TW_TYPE_BOOL8, &cvt()->show_edge_stretch(), "");
			TwAddVarRW(graphics_bar, "- long  rank", TW_TYPE_INT32, &cvt()->stretch_long_rank(), "min=-1 max=100 step=1");
			TwAddVarRW(graphics_bar, "- short rank", TW_TYPE_INT32, &cvt()->stretch_short_rank(), "min=-1 max=100 step=1");
			TwAddVarRW(graphics_bar, "- all   rank", TW_TYPE_INT32, &cvt()->stretch_all_rank(), "min=-1 max=100 step=1");

			//TwAddVarRW(graphics_bar, "show Constrained Edge", TW_TYPE_BOOL8, &cvt()->show_constrained_edge(), "");
			TwAddVarRW(graphics_bar, "show Direction Field", TW_TYPE_BOOL8, &cvt()->show_direction_field(), "");
			///

			TwBar* numerics_bar = TwNewBar("Numerics") ;
            TwDefine(" Numerics position='16 500' size='200 240'") ;  
			TwAddVarRW(numerics_bar, "Density", TW_TYPE_BOOL8, &cvt()->use_density(), "");
			TwAddVarRW(numerics_bar, "nb iter",  TW_TYPE_FLOAT, &nb_iter_, "min=0 max=1000 step=1") ;
            TwAddVarRW(numerics_bar, "Lp order", TW_TYPE_FLOAT, &cvt()->Lp(), "min=0 max=9 step=1") ;
            TwEnumVal aniso_def[] = { {CONSTANT, "const"}, {R_INV, "1/R"}, {BORDER, "Border"} } ;
            TwType tw_aniso = TwDefineEnum("AnisoType", aniso_def, 3) ;
            TwAddVarRW(numerics_bar, "Aniso", tw_aniso, &cvt()->aniso_mode(), "") ;
            TwType tw_center_mode = TwDefineEnum("CenterMode", "centroid, quasi-incenter") ;
            TwAddVarRW(numerics_bar, "CenterMode", tw_center_mode, &cvt()->center_mode(), "") ;
            TwAddVarRW(numerics_bar, "X scale", TW_TYPE_FLOAT, &cvt()->Xscale(), "min=-2 max=2 step=1") ;            
            TwAddVarRW(numerics_bar, "Y scale", TW_TYPE_FLOAT, &cvt()->Yscale(), "min=-2 max=2 step=1") ;            
            TwEnumVal cvt_def[] = { {LLOYD, "Lloyd"}, {NEWTON, "Newton"}, {THETA, "Theta"} } ;
            TwType tw_cvt = TwDefineEnum("CVTType", cvt_def, 3) ;
            TwAddVarRW(numerics_bar, "Mode", tw_cvt, &cvt()->mode(), "") ;
            TwAddVarRW(numerics_bar, "Bndry.", TW_TYPE_BOOL8, &cvt()->insert_boundary(), "") ;            
            TwAddVarRW(numerics_bar, "Edit", TW_TYPE_BOOL8, &edit_, "") ;    

			TwBar* editing_bar = TwNewBar("Editing") ;
			TwDefine(" Editing position='232 10' size='200 240'") ;            
			TwAddVarRW(editing_bar, "Map edit", TW_TYPE_BOOL8, &cvt()->map_edit_mode(), "") ;    
			TwAddButton(editing_bar, "delaunay 2 mesh", CB_convert_mesh, 0, 0);
			TwAddButton(editing_bar, "create grid", CB_create_grid, 0, 0);
			TwAddVarRW(editing_bar, "show mesh.", TW_TYPE_BOOL8, &cvt()->show_mesh(), "") ;            
			TwEnumVal edit_def[] = {{EDGE_FLIP, "edge flip"}, {EDGE_COLLAPSE, "edge collapse"}, {VERTEX_SPLIT, "vertex split"},
			{V4_SPLIT, "v4 split"}, {V4_CREATE, "v4 create"}, {VL7_SPLIT, "vl7 split"}, {VL7_CREATE, "vl7 create"} } ;
            TwType tw_edit = TwDefineEnum("Edit Mode", edit_def, 7) ;
			TwAddVarRW(editing_bar, "Edit Mode", tw_edit, &cvt()->editor()->edit_mode(), "") ;
			TwAddButton(editing_bar, "Action", CB_action, 0, 0);
			
			///dxy add: edit bar
			TwBar* edit_bar = TwNewBar("Edit Topology");
			TwAddVarRW(edit_bar, "Edit", TW_TYPE_BOOL8, &edit_, "") ;  //move to here
			TwAddVarRW(edit_bar, "select", TW_TYPE_BOOL8, &select_, ""); //like edit_
			TwType tw_select_mode = TwDefineEnum("SelectMode", "vertex, edge, face, separatrix");
			TwAddVarRW(edit_bar, "SelectMode", tw_select_mode, &cvt()->delaunay()->select_mode(), "");
			TwAddVarRW(edit_bar, "show selected", TW_TYPE_BOOL8, &cvt()->show_selected(), "");
			TwAddSeparator(edit_bar, "Short_Long", "");
			TwAddVarRW(edit_bar, "Long Edge", TW_TYPE_BOOL8, &cvt()->show_long_edge(), ""); 
			TwAddVarRW(edit_bar, "Short Edge", TW_TYPE_BOOL8, &cvt()->show_short_edge(), ""); 
			TwAddSeparator(edit_bar, "5-7 Pair", "");
			TwAddVarRW(edit_bar, "Num Pair75 Move", TW_TYPE_INT32, &cvt()->num_Pair57_move(), "min=0 max=100 step=1");
			TwAddVarRW(edit_bar, "Num Smooth", TW_TYPE_INT32, &cvt()->num_Pair57_smoothing(), "min=0 max=50 step=1");
			//TwAddSeparator(edit_bar, "Topo Delaunay Optimize", "");
			//TwAddVarRW(edit_bar, "Delaunay Optimize", TW_TYPE_BOOL8, &cvt()->delaunay()->use_topo_delaunay_optimize(), "");
			///

			///dxy add: direction bar
			TwBar* direction_bar = TwNewBar("Direction Field") ;
			TwAddVarRW(direction_bar, "Direction Field", TW_TYPE_BOOL8, &cvt()->direction_field(), "");
			TwType tw_edge_weight_mode = TwDefineEnum("EdgeWeightMode", "dual_edge, lloyd_energy");
			TwAddVarRW(direction_bar, "Edge Weight", tw_edge_weight_mode, &cvt()->direction_edge_weight_mode(), "");
			TwType tw_direction_field_mode = TwDefineEnum("FieldMode", "constant, circle, sin");
			TwAddVarRW(direction_bar, "Field Type", tw_direction_field_mode, &cvt()->direction_field_mode(), "");
			TwAddVarRW(direction_bar, "Direction Angle", TW_TYPE_FLOAT, &cvt()->direction_angle(), "min=-1.0 max=1.0 step=0.0001");
			TwAddVarRW(direction_bar, "Direction Factor", TW_TYPE_FLOAT, &cvt()->direction_factor(), "min=0.0 max=10.0 step=0.01");
			///
			TwAddSeparator(direction_bar, "Quad", "");
			TwAddVarRW(direction_bar, "Quad Field", TW_TYPE_BOOL8, &cvt()->use_quad_field(), "");
			TwAddVarRW(direction_bar, "Quad Factor", TW_TYPE_FLOAT, &cvt()->quad_factor(), "min=0.0 max=1.0 step=0.1");
			///
			TwAddSeparator(direction_bar, "Topo", "");
			TwAddVarRW(direction_bar, "Topo Optimization", TW_TYPE_BOOL8, &cvt()->use_topo_optimization(), "");
			TwAddVarRW(direction_bar, "Balanced", TW_TYPE_BOOL8, &cvt()->use_balanced_stretch(), "");
			//
			TwAddSeparator(direction_bar, "Save", "");
			TwAddVarRW(direction_bar, "Auto Save pts", TW_TYPE_BOOL8, &cvt()->use_auto_save(), "");
			// 
            viewer_properties_->add_separator("Graphics") ;
            viewer_properties_->add_slider("Vertices", cvt()->vertices_size()) ;
            viewer_properties_->add_slider("Centers", cvt()->centers_size()) ;
            viewer_properties_->add_toggle("Primal mesh", cvt()->show_primal_mesh()) ;
            viewer_properties_->add_toggle("Dual mesh", cvt()->show_dual_mesh()) ;
			viewer_properties_->add_toggle("Colorize", cvt()->colorize()) ;
			viewer_properties_->add_toggle("Snap border", cvt()->snap_boundary()) ;
            viewer_properties_->add_toggle("Non-hex", cvt()->show_cells()) ;
            viewer_properties_->add_toggle("Energy", cvt()->show_energy()) ;
            viewer_properties_->add_toggle("Bkgnd. field", cvt()->show_field()) ;
            viewer_properties_->add_slider("Quads", cvt()->quad_ratio(), 0.1, 1.5) ;
			
            viewer_properties_->add_separator("Optimizer") ;
            viewer_properties_->add_slider("Lp order", cvt()->Lp(), 0.0, double(DelaunayCVT::MAX_P))->set_integer(GL_TRUE) ;
            viewer_properties_->add_enum("Aniso", cvt()->aniso_mode(), GlutViewerGUI::LabelList() | AnisoModeNames) ;            
            viewer_properties_->add_slider("X scale", cvt()->Xscale(), -2.0, 2.0)->set_integer(GL_TRUE) ;
            viewer_properties_->add_slider("Y scale", cvt()->Yscale(), -2.0, 2.0)->set_integer(GL_TRUE) ;
            viewer_properties_->add_enum("Mode", cvt()->mode(), GlutViewerGUI::LabelList() | CVTModeNames) ;
            viewer_properties_->add_toggle("Bndry.", cvt()->insert_boundary()) ;
            viewer_properties_->add_toggle("Edit", edit_) ;

            toggle_skybox_CB() ;

            glut_viewer_add_toggle('B', glut_viewer_is_enabled_ptr(GLUT_VIEWER_BACKGROUND), "switch Color/BW") ;
//            glut_viewer_add_toggle('T', &viewer_properties_->visible(), "viewer properties") ;
            viewer_properties_->hide() ;
            glut_viewer_add_toggle('T', glut_viewer_is_enabled_ptr(GLUT_VIEWER_TWEAKBARS), "Tweak bars") ;

            glut_viewer_disable(GLUT_VIEWER_BACKGROUND) ;
        }

        virtual GLboolean mouse(float x, float y, int button, enum GlutViewerEvent event) {
            static int timestamp = 0 ;
            static int last_timestamp = 0 ;
            timestamp++ ;
            static int mode = 0 ;

            GLdouble p[3] ;
            GLdouble v[3] ;
            if(GeexApp::mouse(x, y, button, event)) { return GL_TRUE ; }
            if(edit_) {
 
				bool change = false ;
                if(event == GLUT_VIEWER_DOWN) { mode = button ; change = true ; }
                
                if(event == GLUT_VIEWER_UP) { return GL_FALSE ; }

                glut_viewer_get_picked_ray(p,v) ;
                vec2 pt(p[0], p[1]) ;

				//if(mode == 0 && (change || timestamp - last_timestamp > 3)) {
				//	if(cvt()->in_boundary(pt)) {
				//		bool suc = cvt()->pick_vertex(pt) ;
				//		if(!suc) {
				//			cvt()->pick_edge(pt) ;
				//		}
				//		return GL_TRUE ;
				//	}
				//}

                if(mode == 0 && (change || timestamp - last_timestamp > 3)) {
                    if(cvt()->in_boundary(pt)) {
						if(cvt()->period()) cvt()->clear_copies() ;
                        cvt()->begin_insert() ;
                        cvt()->insert(pt) ;
                        cvt()->end_insert() ;
						if(cvt()->period()) {
							cvt()->insert_copies(cvt()->pvd_mode()==FULL_COPY, false) ;
						}
                    }
                    last_timestamp = timestamp ;
                }

                if(mode == 1 && (change || timestamp - last_timestamp > 3)) {
                    if(cvt()->in_boundary(pt)) {
						if(cvt()->period()) cvt()->clear_copies() ;
                        //cvt()->begin_insert() ;
//                         cvt()->remove(pt) ;
//                         cvt()->end_insert(false) ;
//                         cvt()->begin_insert() ;
//                         cvt()->insert(pt) ; // ->locked = true ;
						cvt()->move_nearest_vertex_to(pt);
                        //cvt()->end_insert() ;
						if(cvt()->period()) cvt()->insert_copies(cvt()->pvd_mode()==FULL_COPY, false) ;
                    }
                    last_timestamp = timestamp ;
                }

                if(mode == 2 && (change || timestamp - last_timestamp > 3)) {
                    if(cvt()->in_boundary(pt)) {
						if(cvt()->period()) cvt()->clear_copies() ;
                        cvt()->begin_insert() ;
                        cvt()->remove(pt) ;
                        cvt()->end_insert() ;
						if(cvt()->period()) cvt()->insert_copies(cvt()->pvd_mode()==FULL_COPY) ;
                    }
                    last_timestamp = timestamp ;
                }

                return GL_TRUE ;
			} 
			// pick vertex, edge
			if(cvt()->map_edit_mode()) {
				bool change = false ;
                if(event == GLUT_VIEWER_DOWN) { mode = button ; change = true ; }
                if(event == GLUT_VIEWER_UP) { return GL_FALSE ; }

                glut_viewer_get_picked_ray(p,v) ;
                vec2 pt(p[0], p[1]) ;

				if(!cvt()->in_boundary(pt)) { return GL_FALSE ; }

				if(mode == 0 && (change || timestamp - last_timestamp > 3)) {
					cvt()->editor()->pick_vertex(pt) ;
					return GL_TRUE ;
				}
				else if(mode == 1 && (change || timestamp - last_timestamp > 3)) {
					cvt()->editor()->pick_edge(pt) ;
					return GL_TRUE ;
				}
				else if(mode == 2 && (change || timestamp - last_timestamp > 3)) {
					cvt()->editor()->clear_picked() ;
					return GL_TRUE ;
				}
			}
			///dxy add: select
			if(!edit_ && select_) {
				bool change = false ;
				if(event == GLUT_VIEWER_DOWN) { mode = button ; change = true ; }
				if(event == GLUT_VIEWER_UP) { return GL_FALSE ; }
				glut_viewer_get_picked_ray(p,v) ;
				vec2 pt(p[0], p[1]) ;
				if(mode == 0 && (change || timestamp - last_timestamp > 3)) { //left click - select/unselect
					if(cvt()->in_boundary(pt)) {
						cvt()->update_select(pt);
					}
					last_timestamp = timestamp ;
				}
				return GL_TRUE ;
			}
			///dxy add end
			return GL_FALSE ;
        }

    private:
        std::string boundary_filename_ ;
        int nb_points_ ;
		int min_max_ ;
        GLfloat nb_iter_ ;
        GLboolean non_convex_ ;
        GLboolean edit_ ;
		///dxy add
		GLboolean select_;
		///
		GLboolean insert_boundary_ ;
		int nb_runs_ ;
    } ;
}

Geex::CVTApp* cvt_app() { return static_cast<Geex::CVTApp*>(Geex::GeexApp::instance()) ; }

void Lloyd() {
    cvt_app()->Lloyd() ;
}

void Lloyd1() {
    cvt_app()->Lloyd(1) ;
}

///dxy add
void Update_Energy_Grad() {
	cvt_app()->Update_Energy_Grad();
}

void Update_Total_Area() {
	cvt_app()->Update_Total_Area();
}

void edge_split() {
	cvt_app()->edge_split();
}

void edge_collapse() {
	cvt_app()->edge_collapse();
}

void edge_flip() {
	//cvt_app()->edge_flip();
	cvt_app()->cvt()->flip_selected_edge_via_map();
}

void pair57_move() {
	cvt_app()->pair57_move();
}

//dxy test
void test() {
	//cvt_app()->cvt()->editor()->test();
	//cvt_app()->cvt()->flip_selected_edge_via_map();
	cvt_app()->cvt()->topo_edit_via_map();
}
///

void NewtonLloyd() {
    cvt_app()->NewtonLloyd() ;
	cvt_app()->cvt()->lloyd_energy() ; // update regularity
}

void reset() {
    cvt_app()->reset() ;
	cvt_app()->cvt()->lloyd_energy() ;
}

void insert_copies() {
	cvt_app()->insert_copies() ;
	cvt_app()->cvt()->lloyd_energy() ;
}

void clear_copies() {
	cvt_app()->clear_copies() ;
}

void euclidean() {
	cvt_app()->cvt()->show_pvd_euclidean() = !cvt_app()->cvt()->show_pvd_euclidean() ;
}
void inc_Lp() {
    cvt_app()->cvt()->Lp() += 1.0 ;
    cvt_app()->cvt()->Lp() = Geex::gx_min(cvt_app()->cvt()->Lp(), float(Geex::DelaunayCVT::MAX_P)) ;
    cvt_app()->cvt()->lp_shader() = int(cvt_app()->cvt()->Lp()) ;
}

void dec_Lp() {
    cvt_app()->cvt()->Lp() -= 1.0 ;
    cvt_app()->cvt()->Lp() = Geex::gx_max(cvt_app()->cvt()->Lp(), 0.0f) ;
    cvt_app()->cvt()->lp_shader() = int(cvt_app()->cvt()->Lp()) ;
}

void save() {
	cvt_app()->save() ;
}

void load() {
	cvt_app()->load() ;
	cvt_app()->cvt()->lloyd_energy() ; // update regularity
}

void stats() { // statistic 
//	cvt_app()->stats() ;
//	cvt_app()->degree_stats(100) ;
	cvt_app()->analysis_pointsets() ;
}

void show_edge_length() {
	cvt_app()->cvt()->compute_edge_length() ;
}

void do_perturb() {
	cvt_app()->cvt()->do_perturb() ;
}

void compute_hessian() {
	cvt_app()->cvt()->compute_hessian() ;
}

void generate_poisson_disk() {
	cvt_app()->generate_poisson_disk() ;
	//std::string fname = "D:/project/gx_pcvt/gx_pcvt2d/pointsets/PVoronoi" ;
	//Geex::psa(true, false, false, fname) ;
}

void update_gapss() {
	cvt_app()->cvt()->update_gaps(cvt_app()->cvt()->sample_radius()) ;
}

void smooth() {
	cvt_app()->cvt()->editor()->smooth() ;
}

void batch_snapshot() {
	std::string dir = Geex::FileSystem::get_project_root() + "/gx_pcvt2d/" ;
	char filename[1024] ;
	int N = 0, nSol = 0 ;
	std::cout << "Please specify N and nSol:" << std::endl;
	std::cout << "N:";
	std::cin >> N;
	std::cout << "nSol:";
	std::cin >> nSol;
	for(int i=0; i<nSol; ++i) {
		sprintf(filename, "%sdata/eg18/%d/%d.pts", dir.c_str(), N, i) ;
		cvt_app()->cvt()->load(filename) ;
		cvt_app()->cvt()->insert_copies(true) ;
		glut_viewer_redraw() ;
		sprintf(filename, "%sdata/eg18/%d/%d.png", dir.c_str(), N, i) ;
		cvt_app()->cvt()->snapshot(filename) ;
	}
}



void TW_CALL CB_action(void *clinetData) {
	cvt_app()->cvt()->editor()->do_editing() ;
}

void TW_CALL CB_convert_mesh(void *clientData) {
	cvt_app()->cvt()->map_edit_mode() = true ;
	cvt_app()->cvt()->editor()->delaunay_to_map() ;
	cvt_app()->cvt()->show_mesh() = true ;
	cvt_app()->cvt()->show_vertices() = false;
	cvt_app()->cvt()->show_primal_mesh() = false ;
	cvt_app()->cvt()->show_dual_mesh() = false ;
}

void TW_CALL CB_create_grid(void *clientData) {
	double r = cvt_app()->cvt()->sample_radius() ;
//	cvt_app()->cvt()->editor()->create_grid_map(r) ;
	cvt_app()->cvt()->editor()->create_regular_map(r) ;
	cvt_app()->cvt()->show_mesh() = true ;
}

void Lloyd_fpo() {
	cvt_app()->Lloyd_fpo() ;
}

void insert_grid() {
	cvt_app()->insert_grid() ;
}

int main(int argc, char** argv) {
    Geex::initialize() ;
    Geex::CVTApp app(argc, argv) ;
    glut_viewer_add_key_func('k', Lloyd, "Lloyd iterations") ;
    glut_viewer_add_key_func('K', Lloyd1, "Lloyd one iteration") ;
	glut_viewer_add_key_func('f', Lloyd_fpo, "Lloyd fpo iterations") ;
    glut_viewer_add_key_func('m', NewtonLloyd, "Newton-Lloyd iterations") ;
	glut_viewer_add_key_func('M', smooth, "smooth point set") ;
    glut_viewer_add_key_func('Z', reset, "reset") ;
	glut_viewer_add_key_func('g', insert_grid, "insert_grid") ;
	glut_viewer_add_key_func('j', compute_hessian, "Hessian") ;
	glut_viewer_add_key_func('i', insert_copies, "Insert copies") ;
	glut_viewer_add_key_func('c', clear_copies, "clear copies") ;
	glut_viewer_add_key_func('u', euclidean, "show euclidean") ;
	glut_viewer_add_key_func('s', save, "save points") ;
	glut_viewer_add_key_func('o', load, "load points") ;
	glut_viewer_add_key_func('d', stats, "statistic") ;
	glut_viewer_add_key_func('z', generate_poisson_disk, "generate poisson disk") ;
	glut_viewer_add_key_func('v', update_gapss, "update void regions") ;
	glut_viewer_add_key_func('b', batch_snapshot, "snapshots") ;
	glut_viewer_add_key_func('p', do_perturb, "perturb") ;
    glut_viewer_add_toggle('e', &(cvt_app()->edit()), "edit") ;
    glut_viewer_add_key_func('<', dec_Lp, "decrement Lp") ;
    glut_viewer_add_key_func('>', inc_Lp, "increment Lp") ;
	///dxy add
	glut_viewer_add_key_func('n', Update_Energy_Grad, "Update energy and grad");
	glut_viewer_add_key_func('A', Update_Total_Area, "Update total area");
	glut_viewer_add_key_func('S', edge_split, "edge split");
	glut_viewer_add_key_func('C', edge_collapse, "edge collapse");
	glut_viewer_add_key_func('F', edge_flip, "edge flip");
	glut_viewer_add_key_func('P', pair57_move, "5-7 pair move");
	//dxy test
	glut_viewer_add_key_func(',', test, "test");
	///
    glut_viewer_disable(GLUT_VIEWER_3D) ;
    app.main_loop() ;
    Geex::terminate() ;
}
