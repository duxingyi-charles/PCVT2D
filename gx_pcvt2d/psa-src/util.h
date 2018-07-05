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

#ifndef UTIL_H
#define UTIL_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cfloat>

#define PI_2    1.57079632679489661923f
#define PI      3.14159265358979323846f
#define TWOPI   6.28318548202514648437f
#define SQRT1_2 0.70710678118654752440f
#define SQRT2   1.41421356237309504880f
#define SQRT3   1.73205080756887729352f


/**
 *  Clamps a value to a minimum of 0.0
 */
inline void clamp0(float &f) {
    f = (f < 0.0f ? 0.0f : f);
}

/**
 *  Clamps a value to a maximum of 1.0
 */
inline void clamp1(float &f) {
    f = (f > 1.0f ? 1.0f : f);
}

/**
 *  Clamps a value to be in the interval [0, 1]
 */
inline void clamp01(float &f) {
    clamp0(f);
    clamp1(f);
}

/**
 *  Squares the specified value
 */
inline float sqr(float f) {
    return f * f;
}

/**
 *  Returns the least power of 2 greater than or equal to x
 */
inline unsigned nextpow2(unsigned x) {
    x = x - 1;
    x = x | (x >> 1);
    x = x | (x >> 2);
    x = x | (x >> 4);
    x = x | (x >> 8);
    x = x | (x >> 16);
    return x + 1;
}

/**
 *    Takes the logarithm to base 2
 */
inline float log2f(float f) {
    return log10f(f) / log10f(2.0f);
}

/**
 *  Returns f in decibel
 */
inline float decibel(float f) {
    return 10.0f * log10f(f);
}

/**
 *  Sets all values of the specified buffer to 0
 */
inline void SetZero(float *&buf, int size) {
    for (int i = 0; i < size; ++i) {
        buf[i] = 0.0f;
    }
}

/**
 *  Divide all values of the specified buffer by the specified value
 */
inline void Divide(float *&buf, int size, float d) {
    assert(d != 0.0f);
    for (int i = 0; i < size; ++i) {
        buf[i] /= d;
    }
}

/**
 *  Converts all values of the specified buffer to decibel
 */
inline void Decibel(float *&buf, int size) {
    for (int i = 0; i < size; ++i) {
        buf[i] = decibel(buf[i]);
    }
}

/**
 *  Determines if two values are equivalent considering a small interval
 */
#define EQUIV_EPSILON    (1e-9)
inline bool equivalent(float f1, float f2) {
    return fabsf(f1 - f2) < EQUIV_EPSILON;
}

/**
 *  Convert float from little endian to big endian
 */
inline float big_endian(float f) {
    union {
        float f;
        char c[4];
    } a, b;
    a.f = f;
    b.c[0] = a.c[3];
    b.c[1] = a.c[2];
    b.c[2] = a.c[1];
    b.c[3] = a.c[0];
    return b.f;
}

/**
 *  Returns the squared euclidean distance between two points given by
 *  (u, v) and (x, y) considering the unit square.
 */
inline float dist_sqr(float u, float v, float x, float y) {
    return sqr(u - x) + sqr(v - y);
}

/**
 *  Returns the squared euclidean distance between two points given by
 *  (u, v) and (x, y) considering the unit torus.
 */
inline float dist_sqr_torus(float u, float v, float x, float y) {
    float minX, minY;
    minX = fabsf(u - x);
    minY = fabsf(v - y);
    minX = (minX > 0.5f) ? 1.0f - minX : minX;
    minY = (minY > 0.5f) ? 1.0f - minY : minY;
    return sqr(minX) + sqr(minY);
}

/**
 *  Collects basic statistical data
 */
struct stats {
    float mean;
    float sd;
    float min;
    float max;
};

/**
 *  Collects information about the plot formatting
 */
struct info {
    double color[4];    // Graph RGBA color
    int start;          // Graph starting value (useful for excluding DC peak)
    double yscale[2];   // Min and max for y-axis
    double ytrange[2];  // Tick range for y-axis
    double ytsteps[2];  // Small and large tick spacing for y-axis
    double ref;         // Value for reference axis
    double nlvl;        // Value for background noise (0.0 equals no noise)
};

const info info_rp = {  // Default radial power formatting
    { 0.0, 0.0, 0.0, 1.0 },
    1,
    { -0.1, 2.1 },
    { -0.1, 2.1 },
    {  0.1, 0.5 },
    1.0,
    0.0
};
const info info_ani = {  // Default anisotropy formatting
    { 0.0, 0.0, 0.0, 1.0 },
    4,
    { -12.5, 12.5 },
    { -12.0, 12.0 },
    { 1.0, 5.0 },
    0.0,
    -10.0
};

/**
 *  Computes anisotropy formatting parameters for a given number of point sets
 */
inline info ani_format(int nsets) {
    info i = info_ani;
    const double d = decibel(nsets);
    i.nlvl = -d;
    if (nsets > 10) {
        i.yscale[0] = -1.25 * d;
        i.yscale[1] =  1.25 * d;
        i.ytrange[0] = -floor(1.2 * d);
        i.ytrange[1] =  floor(1.2 * d);
        i.ytsteps[0] =  floor(d / 10.0);
        i.ytsteps[1] =  floor(d /  2.0);
    }
    return i;
}

#endif    // UTIL_H

