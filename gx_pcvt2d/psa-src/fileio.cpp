/**
 * This file is part of the point set analysis tool psa
 * 
 * Copyright 2011, Thomas Schl\"{o}mer, thomas.schloemer@uni-konstanz.de
 * 
 * psa is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "fileio.h"

#include <iostream>
#include <sstream>

#include <cairo/cairo.h>
#include <cairo/cairo-pdf.h>
#include <cairo/cairo-ps.h>

#include "util.h"

using namespace std;


enum text_align {
    text_align_left,
    text_align_center,
    text_align_right
};


void FileList(string fname, list<string> &fnames)
{
    fnames.clear();
    
    if ( FolderExists(fname) )
    {
        string dname = string(fname);
        RemoveTrailingSlash(dname);
        
        DirIterOS dir(dname);
        DirEntry entry;
        
        while ( dir.get_next(entry) )
            if ( !entry.is_dir )
                if ( (entry.name.find(FILE_EXT_RPS) != string::npos) ||
                     (entry.name.find(FILE_EXT_TXT) != string::npos) )
                    fnames.push_back(dname + "/" + entry.name);
    } else {
        if ( (fname.find(FILE_EXT_RPS) != string::npos) ||
             (fname.find(FILE_EXT_TXT) != string::npos) ) {
            fnames.push_back(fname);
        }
    }
}


bool ReadNumPointsRPS(string fname, int &npoints)
{
    if ( fname.find(FILE_EXT_RPS) == string::npos ) {
        cerr << "ReadNumPointsRPS: No " << FILE_EXT_RPS << " file " << fname
        << endl;
        return false;
    }
    
    FILE *fp;
    if ( !(fp = fopen(fname.c_str(), "rb")) ) {
        cerr << "ReadNumPointsRPS: Error opening " << fname << endl;
        return false;
    }
    
    long filesize;    
    fseek(fp, 0, SEEK_END); 
    filesize = ftell(fp); 
    fseek(fp, 0, SEEK_SET);
    
    npoints = filesize / (2 * sizeof(float));
    
    fclose(fp);
    return true;
}


bool ReadNumPointsTXT(string fname, int &npoints)
{
    if ( fname.find(FILE_EXT_TXT) == string::npos ) {
        cerr << "ReadNumPointsTXT: No " << FILE_EXT_TXT << " file "
             << fname << endl;
        return false;
    }
    
    FILE *fp;
    if ( !(fp = fopen(fname.c_str(), "r")) ) {
        cerr << "ReadNumPointsTXT: Error opening " << fname << endl;
        return false;
    }
    
    fscanf(fp, "%d\n", &npoints);

    fclose(fp);
    return true;
}


bool ReadNumPoints(string fname, int &npoints)
{
    if ( fname.find(FILE_EXT_RPS) != string::npos ) {
        return ReadNumPointsRPS(fname, npoints);
    }
    
    if ( fname.find(FILE_EXT_TXT) != string::npos ) {
        return ReadNumPointsTXT(fname, npoints);
    }
    
    cerr << "No " << FILE_EXT_RPS << " or " << FILE_EXT_TXT << " file "
         << fname << endl;
    return false;
}


bool LoadPointsRPS(string fname, float **points, int &npoints)
{
    if ( !FileExists(fname) ) {
        cerr << "LoadRPS: File does not exist " << fname << endl;
        return false;
    }
    
    if ( (fname.find(FILE_EXT_RPS) == string::npos) ) {
        cerr << "LoadRPS: No " << FILE_EXT_RPS << " file " << fname << endl;
        return false;
    }
    
    FILE *fp;
    if ( !(fp = fopen(fname.c_str(), "r")) ) {
        cerr << "LoadRPS: Error opening " << fname << endl;
        return false;
    }
    
    long filesize;    
    fseek(fp, 0, SEEK_END); 
    filesize = ftell(fp); 
    fseek(fp, 0, SEEK_SET);
    
    npoints = filesize / (2 * sizeof(float));
    *points = new float[npoints * 2];
    
    if ( !fread(*points, sizeof(float), npoints * 2, fp) ) {
        cerr << "LoadRPS: Cannot read point set from " << fname << endl;
        fclose(fp);
        return false;
    }
    
    fclose(fp);
    return true;
}


bool LoadPointsTXT(string fname, float **points, int &npoints)
{
    if ( fname.find(FILE_EXT_TXT) == string::npos ) {
        cerr << "LoadTXT: No " << FILE_EXT_TXT << " file " << fname << endl;
        return false;
    }
    
    FILE* fp;
    if ( !(fp = fopen(fname.c_str(), "r")) ) {
        cerr << "LoadTXT: Error opening " << fname << endl;
        return false;
    }
    
    fscanf(fp, "%d\n", &npoints);
    *points = new float[npoints * 2];
    
    for (int i = 0; i < npoints; ++i) {
        fscanf(fp, "%f %f\n", &(*points)[2*i], &(*points)[2*i+1]);
    }
    
    fclose(fp);
    return true;
}


bool LoadPoints(string fname, float **points, int &npoints)
{
    if ( fname.find(FILE_EXT_RPS) != string::npos ) {
        return LoadPointsRPS(fname, points, npoints);
    }
    
    if ( fname.find(FILE_EXT_TXT) != string::npos ) {
        return LoadPointsTXT(fname, points, npoints);
    }
    
    cerr << "No " << FILE_EXT_RPS << " or " << FILE_EXT_TXT << " file "
         << fname << endl;
    return false;
}


bool SaveRPS(string fname, const float * const points, int npoints)
{
    FILE *fp;
    if ( !(fp = fopen(fname.c_str(), "wb")) ) {
        cerr << "SaveRPS: Cannot write point set to " << fname << endl;
        return false;
    }
    
    if ( !fwrite(points, sizeof(float), npoints * 2, fp) ) {
        fclose(fp);
        return false;
    }
    
    cout << npoints << " points written to " << fname << endl;
    
    fclose(fp);
    return true;
}


bool SaveTXT(string fname, const float * const points, int npoints)
{
    FILE *fp;
    if ( !(fp = fopen(fname.c_str(), "w")) ) {
        cerr << "SaveTXT: Cannot write point set to " << fname << endl;
        return false;
    }
    
    fprintf(fp, "%d\n", npoints);
    
    for (int i = 0; i < npoints; ++i) {
        fprintf(fp, "%f %f\n", points[2*i], points[2*i+1]);
    }
    
    fclose(fp);
    return true;
}


bool SaveFloatPNG(string fname, const float * const image, int resX, int resY,
                  bool flipped)
{
    int stride = cairo_format_stride_for_width(CAIRO_FORMAT_RGB24, resX);
    unsigned char *data = new unsigned char[stride * resY];
    
    for (int x = 0; x < resX; ++x)
    {
        for (int y = 0; y < resY; ++y)
        {
            int flipped_y = (flipped ? resY - 1 - y : y);
            
            float f = image[x + flipped_y*resX];
            clamp01(f);
            
            data[4*(x + flipped_y*resX)  ] = (unsigned char)(f * 255.0f);
            data[4*(x + flipped_y*resX)+1] = (unsigned char)(f * 255.0f);
            data[4*(x + flipped_y*resX)+2] = (unsigned char)(f * 255.0f);
            // Fourth component stays empty
        }
    }

    cairo_surface_t *surface = cairo_image_surface_create_for_data(data,
                                       CAIRO_FORMAT_RGB24, resX, resY,
                                       stride);
    cairo_status_t status = cairo_surface_write_to_png(surface,
                                                       fname.c_str());
    cairo_surface_destroy(surface);
    delete []data;

    return (status == CAIRO_STATUS_SUCCESS);
}


void PrintToContext(cairo_t *cr, string label, double atX, double atY,
                    string family, text_align align)
{
    cairo_text_extents_t extents;
    cairo_select_font_face(cr, family.c_str(), CAIRO_FONT_SLANT_NORMAL,
                           CAIRO_FONT_WEIGHT_NORMAL);
    cairo_text_extents(cr, label.c_str(), &extents);
    
    double x, y;
    
    if (align == text_align_left)
        x = atX - extents.x_bearing;
    else if (align == text_align_center)
        x = atX - (extents.width / 2.0 + extents.x_bearing);
    else // (align == text_align_right)
        x = atX - (extents.width + extents.x_bearing);
    
    y = atY - (extents.height/2 + extents.y_bearing);
    
    cairo_move_to(cr, x, y);
    cairo_show_text(cr, label.c_str());
    cairo_stroke(cr);
}


void DrawPointsToContext(cairo_t *cr, double atX, double atY,
                         const float * const points, int npoints,
                         int width, int height)
{
    const float radius = 2.0;
    cairo_set_source_rgba(cr, 0.0, 0.0, 0.0, 1.0);
    
    for (int i = 0; i < npoints; ++i) {
        double x = points[2*i] * width;
        double y = (1.0 - points[2*i+1]) * height;
        cairo_arc(cr, atX + x, atY + y, radius, 0, PI * 2.0);
        cairo_fill(cr);
    }
}


void SavePointPlot(string fname, bool eps, const float * const points,
                   int npoints, int width, int height)
{
    cairo_surface_t *surface;
    
    if (eps) {
        surface = cairo_ps_surface_create(fname.c_str(), width, height);
        cairo_ps_surface_set_eps(surface, true);
    } else
        surface = cairo_pdf_surface_create(fname.c_str(), width, height);
    
    cairo_t *cr = cairo_create(surface);
    
    DrawPointsToContext(cr, 0.0, 0.0, points, npoints, width, height);
    
    cairo_show_page(cr);
    cairo_surface_destroy(surface);
    cairo_destroy(cr);
}


void DrawPlotToContext(cairo_t *cr, float atX, float atY,
                       const float * const buf, size_t size, int width,
                       int height, info i)
{
    double refAxis   = (i.ref  - i.yscale[0]) / (i.yscale[1] - i.yscale[0]);
    double noiseAxis = (i.nlvl - i.yscale[0]) / (i.yscale[1] - i.yscale[0]);
    
    // Draw reference axes
    cairo_set_line_width(cr, 2.0);
    cairo_set_source_rgba(cr, 0.6, 0.6, 0.6, 1.0);
    double dashes[] = { 6.0, 3.0 };
    cairo_set_dash(cr, dashes, 2, 0.0);
    cairo_move_to(cr, atX, atY + (1.0 - refAxis) * height);
    cairo_line_to(cr, atX + width, atY + (1.0 - refAxis) * height);
    cairo_stroke(cr);
    
    if (i.nlvl != 0.0) {
        cairo_move_to(cr, atX, atY + (1.0 - noiseAxis) * height);
        cairo_line_to(cr, atX + width, atY + (1.0 - noiseAxis) * height);
        cairo_stroke(cr);
    }
    
    // Draw plot
    cairo_set_line_width(cr, 1.0);
    cairo_set_dash(cr, NULL, 0, 0.0);
    cairo_set_source_rgba(cr, i.color[0], i.color[1], i.color[2], i.color[3]);
    
    float x, y;
    
    for (int j = i.start; j < (int) size; ++j) {
        x = j / (float) size;
        y = (buf[j] - i.yscale[0]) / (i.yscale[1] - i.yscale[0]);
        clamp01(y);
        y = 1.0 - y;
        
        if (j == i.start)
            cairo_move_to(cr, atX + x * width, atY + y * height);
        else
            cairo_line_to(cr, atX + x * width, atY + y * height);
    }
    cairo_stroke(cr);
}


void SaveSpectralPlot(string fname, bool eps, string label,
                      const float * const buf, size_t size, int width,
                      int height, info i)
{
    double ratio = width / (double) height;
    double canvas_width  = width  * (1.0 + 0.25 / ratio);
    double canvas_height = height * 1.25;
    
    double offset_x = canvas_width - (width * (1.0 + 0.05 / ratio));
    double offset_y = height * 0.05;
    
    double label_fontsize = 28.0;
    double axes_fontsize  = 18.0;
    
    double tick_len_large = 0.015 * std::max(width, height);
    double tick_len_small = tick_len_large / 2.0;

    cairo_surface_t *surface;
    
    if (eps) {
        surface = cairo_ps_surface_create(fname.c_str(), canvas_width,
                                          canvas_height);
        cairo_ps_surface_set_eps(surface, true);
    } else
        surface = cairo_pdf_surface_create(fname.c_str(), canvas_width,
                                           canvas_height);
    
    cairo_t *cr = cairo_create(surface);
    
    // Plot
    DrawPlotToContext(cr, offset_x, offset_y, buf, size, width, height, i);
    
    // Draw frame
    cairo_set_source_rgba(cr, 0.0, 0.0, 0.0, 1.0);
    cairo_set_line_width(cr, 1.0);
    cairo_rectangle(cr, offset_x, offset_y, width, height);
    cairo_stroke(cr);
    
    
    // x-axis captions
    
    // Label
    cairo_set_font_size(cr, label_fontsize);
    PrintToContext(cr, "frequency", offset_x + 0.5 * width,
                   offset_y + height + 0.6 * (canvas_height - height),
                   "Sans", text_align_center);
    cairo_set_font_size(cr, axes_fontsize);
    
    // Ticks
    int nticks = 20;
    int xsteps[2];
    xsteps[0] = size / nticks;
    xsteps[1] = 5 * xsteps[0];

    PrintToContext(cr, "0", offset_x, offset_y + height + 13,
                   "Sans", text_align_center);
    for (int x = xsteps[0]; x < (int) size; x += xsteps[0])
    {
        double x_f = x / (double) size;
        bool isLargeTick = (x % xsteps[1] == 0);
        double tick_len = (isLargeTick ? tick_len_large : tick_len_small);
        
        cairo_move_to(cr, offset_x + x_f * width, offset_y + height);
        cairo_line_to(cr, offset_x + x_f * width,
                      offset_y + height - tick_len);
        cairo_stroke(cr);
        
        if (isLargeTick) {
            stringstream tick_label;
            tick_label << x;
            PrintToContext(cr, tick_label.str(), offset_x + x_f * width,
                           offset_y + height + 13, "Sans", text_align_center);
        }
    }
    
    // Mark critical frequency
    double critical_f = SQRT1_2;
    double dashes[] = { 6.0, 3.0 };
    cairo_set_dash(cr, dashes, 2, 0.0);
    cairo_set_source_rgba(cr, 0.6, 0.6, 0.6, 1.0);
    cairo_move_to(cr, offset_x + critical_f * width, offset_y + height);
    cairo_line_to(cr, offset_x + critical_f * width, offset_y);
    cairo_stroke(cr);
    cairo_set_source_rgba(cr, 0.0, 0.0, 0.0, 1.0);
    cairo_set_dash(cr, NULL, 0, 0.0);
    
    
    // y-axis captions
    
    // Label
    cairo_matrix_t m;
    cairo_matrix_init_identity(&m);
    cairo_matrix_scale(&m, label_fontsize, label_fontsize);
    cairo_matrix_rotate(&m, -PI_2);
    cairo_set_font_matrix(cr, &m);
    PrintToContext(cr, label, 0.25 * offset_x, offset_y + height / 2.0,
                   "Sans", text_align_center);
    
    // Ticks
    cairo_set_font_size(cr, axes_fontsize);
    for (double y = i.ytrange[0]; y <= i.ytrange[1]; y += i.ytsteps[0])
    {
        double y_f = ( y - i.yscale[0] ) / ( i.yscale[1] - i.yscale[0] );
        bool isLargeTick = equivalent( fmodf(y, i.ytsteps[1]), 0.0f );
        double tick_len = (isLargeTick ? tick_len_large : tick_len_small);
        
        cairo_move_to(cr, offset_x, offset_y + (1.0 - y_f) * height);
        cairo_line_to(cr, offset_x + tick_len,
                      offset_y + (1.0 - y_f) * height);
        cairo_stroke(cr);
        
        if (isLargeTick) {
            stringstream tick_label;
            tick_label << y;
            PrintToContext(cr, tick_label.str(), 0.9 * offset_x,
                           offset_y + (1.0 - y_f) * height, "Sans",
                           text_align_right);
        }
    }
    
    cairo_show_page(cr);
    cairo_surface_destroy(surface);
    cairo_destroy(cr);
}


void SaveComposite(std::string fname, bool eps, const float * const points,
                   int npoints, const float * const image, int image_resX,
                   int image_resY, const float * const power,
                   const float * const anisotropy, int nannuli, stats l2disc,
                   stats mindist, stats mindist_norm, stats avgmd_norm,
                   int nsets, int cell_w, int cell_h)
{
    int comp_w = cell_w * 2;
    int comp_h = cell_h * 1.5;
    
    cairo_surface_t* surface;
    
    if (eps) {
        surface = cairo_ps_surface_create(fname.c_str(), comp_w, comp_h);
        cairo_ps_surface_set_eps(surface, true);
    } else
        surface = cairo_pdf_surface_create(fname.c_str(), comp_w, comp_h);
    
    cairo_t* cr = cairo_create(surface);
    
    // Draw points
    DrawPointsToContext(cr, 0.0, 0.0, points, npoints, cell_w, cell_h);
    
    // Draw radial power and anisotropy
    DrawPlotToContext(cr,    0.0, cell_h,      power, nannuli, cell_w,
                      cell_h / 2, info_rp);
    DrawPlotToContext(cr, cell_w, cell_h, anisotropy, nannuli, cell_w,
                      cell_h / 2, ani_format(nsets));
    
    // Mark critical frequencies
    float critical_f = SQRT1_2;
    double dashes[] = { 6.0, 3.0 };
    cairo_set_dash(cr, dashes, 2, 0.0);
    cairo_set_source_rgba(cr, 0.6, 0.6, 0.6, 1.0);
    cairo_move_to(cr, critical_f * cell_w, comp_h);
    cairo_line_to(cr, critical_f * cell_w, comp_h - cell_h / 2);
    cairo_stroke(cr);
    cairo_move_to(cr, cell_w + critical_f * cell_w, comp_h);
    cairo_line_to(cr, cell_w + critical_f * cell_w, comp_h - cell_h / 2);
    cairo_stroke(cr);
    cairo_set_source_rgba(cr, 0.0, 0.0, 0.0, 1.0);
    cairo_set_dash(cr, NULL, 0, 0.0);

    
    // Draw spectrum image in cairo RGB24 format
    int stride = cairo_format_stride_for_width(CAIRO_FORMAT_RGB24,
                                               image_resX);
    unsigned char *data = new unsigned char[stride * image_resY];
    
    for (int x = 0; x < image_resX; ++x) {
        for (int y = 0; y < image_resY; ++y) {
            int flipped_y = image_resY - 1 - y;
            float f = image[x + flipped_y*image_resX];
            clamp01(f);
            data[4*(x + y*image_resX)  ] = (unsigned char)(f * 255.0f);
            data[4*(x + y*image_resX)+1] = (unsigned char)(f * 255.0f);
            data[4*(x + y*image_resX)+2] = (unsigned char)(f * 255.0f);
            // Fourth component stays empty
        }
    }
    cairo_surface_t* image_surface = cairo_image_surface_create_for_data(data,
                          CAIRO_FORMAT_RGB24, image_resX, image_resY, stride);
    
    cairo_translate(cr, cell_w, 0);
    cairo_scale(cr, cell_w / (float) image_resX, cell_h / (float) image_resY);
    cairo_set_source_surface(cr, image_surface, 0, 0);
    cairo_paint(cr);
    
    // Draw separators
    cairo_identity_matrix(cr);
    cairo_set_line_width(cr, 1.0);
    cairo_set_source_rgba(cr, 0.0, 0.0, 0.0, 1.0);
    
    cairo_move_to(cr, 0.0, cell_h);
    cairo_line_to(cr, comp_w, cell_h);
    cairo_stroke(cr);
    
    cairo_move_to(cr, cell_w, 0.0);
    cairo_line_to(cr, cell_w, comp_h);
    cairo_stroke(cr);
    
    // Draw stats box
    bool multiple = (nsets > 1);
    double bsize[] = { (multiple ? 0.95 : 0.34) * cell_w, 0.1275 * cell_h };
    double banchor = 0.0125 * cell_w;
    cairo_set_source_rgba(cr, 0.0, 0.0, 0.0, 0.7);
    cairo_rectangle(cr, banchor, banchor, bsize[0], bsize[1]);
    cairo_fill(cr);
    
    cairo_set_source_rgba(cr, 1.0, 1.0, 1.0, 1.0);
    cairo_set_font_size(cr, 12.0);
    char stat[256];
    if (!multiple) {
        sprintf(stat, "L2-norm discr: %.6f", l2disc.mean);
        PrintToContext(cr, stat, 1.75 * banchor, 0.03 * cell_h, "monospace",
                       text_align_left);
        sprintf(stat, "Mindist      : %.6f", mindist.mean);
        PrintToContext(cr, stat, 1.75 * banchor, 0.06 * cell_h, "monospace",
                       text_align_left);
        sprintf(stat, "Norm. mindist: %.6f", mindist_norm.mean);
        PrintToContext(cr, stat, 1.75 * banchor, 0.09 * cell_h, "monospace",
                       text_align_left);
        sprintf(stat, "Norm. avg. md: %.6f", avgmd_norm.mean);
        PrintToContext(cr, stat, 1.75 * banchor, 0.12 * cell_h, "monospace",
                       text_align_left);
    } else {
        sprintf(stat, "L2-norm discr: mean=%.6f sd=%.6f min=%.6f max=%.6f",
                l2disc.mean, l2disc.sd, l2disc.min, l2disc.max);
        PrintToContext(cr, stat, 1.75 * banchor, 0.03 * cell_h, "monospace",
                       text_align_left);
        sprintf(stat, "Mindist      : mean=%.6f sd=%.6f min=%.6f max=%.6f",
                mindist.mean, mindist.sd, mindist.min, mindist.max);
        PrintToContext(cr, stat, 1.75 * banchor, 0.06 * cell_h, "monospace",
                       text_align_left);
        sprintf(stat, "Norm. mindist: mean=%.6f sd=%.6f min=%.6f max=%.6f",
                mindist_norm.mean, mindist_norm.sd, mindist_norm.min,
                mindist_norm.max);
        PrintToContext(cr, stat, 1.75 * banchor, 0.09 * cell_h, "monospace",
                       text_align_left);
        sprintf(stat, "Norm. avg. md: mean=%.6f sd=%.6f min=%.6f max=%.6f",
                avgmd_norm.mean, avgmd_norm.sd, avgmd_norm.min,
                avgmd_norm.max);
        PrintToContext(cr, stat, 1.75 * banchor, 0.12 * cell_h, "monospace",
                       text_align_left);
    }
    
    // Compilation complete
    cairo_show_page(cr);
    cairo_surface_destroy(image_surface);
    delete []data;
    cairo_surface_destroy(surface);
    cairo_destroy(cr);
}

