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


#include "cvt.h"

namespace Geex {


    static GLint create_contouring_texture() {
        // On the quadro, seems that GLU does not want
        // more than 4096 in tex dim (did not get the 
        // problem on the GeForce, that's wierd...)
        const int factor = 1 ;
        const int tex_size = 4096 * factor ;
        const int width = factor ;
        const int zero_width = 5 * factor * width + 1 ;

        // OpenGL texture upload only works with floats
        // with my driver (therefore we do not use default
        // vec3 type that may use doubles).
        typedef vec3g<float> col ;

        Array1d<col> buffer(2*tex_size) ;
        buffer.set_all(col(1.0f,1.0f,1.0f)) ;

        for(unsigned int i=0; i<tex_size; i++) {
            float s = float(i) / float(tex_size - 1) ;
            buffer[i] = s * col(0.0f,0.0f,0.5f) + (1.0f - s) * col(0.5f, 1.0f, 0.0f) ;
        }

        if(true) {
            for(unsigned int i=0; i<tex_size; i++) {
                if((((i / (2*width))*2) % (128*factor)) == 0) {
                    buffer[i] = col(5.0f, 5.0f, 5.0f) ;
                } 
            }
        }

        GLuint result ;
        glGenTextures(1, &result) ;
        glBindTexture(GL_TEXTURE_2D, result) ;
        gluBuild2DMipmaps(
            GL_TEXTURE_2D, GL_RGB_FLOAT16_ATI, tex_size, 1, 
            GL_RGB, GL_FLOAT, buffer.data()
        ) ;
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
        glBindTexture(GL_TEXTURE_2D, 0) ;
        return result ;
    }

    
    CVT::CVT() : DelaunayCVT(this), DelaunayGraphics(this, this), DelaunayIO(this)  {
        nb_frames_ = 10 ;
        lp_shader_ = 1 ;
        colormap_texture_ = 0 ;
        show_lp_ = false ;
        set_shader("CVT2d") ;
    }

    void CVT::do_draw() {
        DelaunayGraphics::draw() ;
        if(show_lp_) {
            glUseProgramObjectARB(shader_) ;
            glBindTexture(GL_TEXTURE_2D, colormap_texture_) ;
            glUniform1iARB(colormap_id_, 0) ;
            glUniform1fARB(lp_id_, gx_max(1.0f, float(lp_shader_))) ;
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL) ;
            FOR_EACH_VERTEX_DT(Delaunay, this, v) {
                glUniform2fARB(center_id_, v->point().x(), v->point().y()) ;
                vec2 U ; vec2 V ;
                query_anisotropy(vec2(v->point().x(), v->point().y()), U, V) ;
                glUniform2fARB(U_id_, U.x, U.y) ;
                glUniform2fARB(V_id_, V.x, V.y) ;
                DelaunayGraphics::draw_dual_facet(v) ;
            }
            glBindTexture(GL_TEXTURE_2D, 0) ;
        }
    }

    void CVT::pre_draw() {
        Geexob::pre_draw() ;
        if(colormap_texture_ == 0) {
            colormap_texture_ = create_contouring_texture() ;
            colormap_id_ = glGetUniformLocationARB(shader_, "colormap") ;
            center_id_ = glGetUniformLocationARB(shader_, "center") ;
            lp_id_ = glGetUniformLocationARB(shader_, "lp") ;
            U_id_ = glGetUniformLocationARB(shader_, "U") ;
            V_id_ = glGetUniformLocationARB(shader_, "V") ;
        }
    }

    void CVT::get_bbox(
            real& x_min, real& y_min, real& z_min,
            real& x_max, real& y_max, real& z_max
    ) {
        Delaunay::get_bbox(x_min, y_min, z_min, x_max, y_max, z_max) ; 
    }


    void CVT::set_frame(int x) {
        if(frame_ != x) {
            if(mode() == LLOYD) {
                lloyd(3, false) ;
            } else {
                newton_lloyd(3, false) ;
            }
        }
        Geexob::set_frame(x) ;
    }
}
