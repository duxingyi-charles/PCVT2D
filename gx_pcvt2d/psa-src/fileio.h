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

#ifndef FILEIO_H
#define FILEIO_H

#include <string>
#include <list>
#include "fileutil.h"
#include "util.h"

// File format .rps ('raw point set'); raw single-precision floating points
const std::string FILE_EXT_RPS = ".rps";

// Plain text files where the first line contains the number of points
const std::string FILE_EXT_TXT = ".txt";


/**
 *  Generates a list of files from the specified file or folder name.
 */
void FileList(std::string fname, std::list<std::string> &fnames);


/**
 *  Reads the number of points from the specified file. The function
 *  determines the file ending and calls the appropriate reading function.
 */
bool ReadNumPoints(std::string fname, int &npoints);


/**
 *  Loads points from the specified file into the specified points buffer.
 *  Points can currently be stored in two formats: as raw binary data in an
 *  .rps file or in a plain .txt file. RPS files are expected to contain raw
 *  single-precision floating point data. TXT files are expected to follow
 *  the simple format
 *
 *  npoints
 *  point[0].x point[0].y
 *  point[1].x point[1].y
 *  point[2].x point[2].y
 *  ...
 *
 *  This function determines the ending of the specified file and calls the
 *  appropriate loading function for each type.
 */
bool LoadPoints(std::string fname, float **points, int &npoints);


/**
 *  Saves the specified points into the specified file in raw binary format.
 *
 *  You may copy+paste this function into your own project to store your
 *  points in the .rps format.
 */
bool SaveRPS(std::string fname, const float * const points, int npoints);


/**
 *  Saves the specified points into the specified plain text file with the
 *  format
 *
 *  npoints
 *  point[0].x point[0].y
 *  point[1].x point[1].y
 *  point[2].x point[2].y
 *  ...
 *
 *  You may copy+paste this function into your own project to store your
 *  points in this format.
 */
bool SaveTXT(std::string fname, const float * const points, int npoints);


/**
 *  Saves a grayscale PNG image from the specified float data via cairo.
 *
 *  It is assumed that 'image' only represents a single-channel, grayscale
 *  image. If 'flipped' equals 'true', the image is flipped horizontally.
 */
bool SaveFloatPNG(std::string fname, const float * const image, int resX,
                  int resY, bool flipped = true);


/**
 *  Saves the point set to the specified file.
 *
 *  It is assumed that points are in the interval [0.0, 1.0)^2.
 */
void SavePointPlot(std::string fname, bool eps, const float * const points,
                   int npoints, int width, int height);


/**
 *  Saves radial power or anisotropy to the specified file.
 *
 *  Plot and axes are conveniently labeled for a 'paper ready' output.
 *  Formatting specifica are controlled by the struct info for
 *  which there are defaults as defined in the util header.
 */
void SaveSpectralPlot(std::string fname, bool eps, std::string label,
                      const float * const buf, size_t size, int width,
                      int height, info i);


/**
 *  Saves a composite of the point set, the spectrum image, radial
 *  power and anisotropy to the specified file in PDF format.
 *
 *  Plots and axes are not labeled but marked with the critical
 *  frequency and reference axes for a quick overview. The complete
 *  composite will have a size of 2*cell_w x 1.5*cell_h.
 */
void SaveComposite(std::string fname, bool eps, const float * const points,
                   int npoints, const float * const image, int image_resX,
                   int image_resY, const float * const power,
                   const float * const anisotropy, int nannuli, stats l2disc,
                   stats mindist, stats mindist_norm, stats avgmd_norm,
                   int nsets, int cell_w, int cell_h);


#endif // FILEIO_H

