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

#include "spectralanalysis.h"

#include <fftw3.h>
#include "util.h"

#ifdef _OPENMP
#include <omp.h>
#endif


void CFT(const float * const points, int npoints, float *&ft, int domain)
{
    const int domain2 = domain / 2;
    
#ifdef _OPENMP
#pragma omp parallel
#endif
{
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (int y = 0; y < domain; ++y) {
        for (int x = 0; x < domain; ++x) {
            float f_x = 0.0f, f_y = 0.0f;
            float w_x = x - domain2;
            float w_y = y - domain2;
            for (int i = 0; i < npoints; ++i) {
                float exp = -TWOPI * (w_x * points[2*i] +
                                      w_y * points[2*i+1]);
                f_x += cosf(exp);
                f_y += sinf(exp);
            }
            ft[2*(x + y*domain)  ] = f_x;
            ft[2*(x + y*domain)+1] = f_y;
        }
    }
}
}


void FFT(const float * const points, int npoints, float *&ft, int domain)
{
    const int dftdomain = std::max(4 * domain, 512);
    
    fftw_complex *in;
    fftw_complex *out;
    
    in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * sqr(dftdomain));
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * sqr(dftdomain));
    
    for (int c = 0; c < sqr(dftdomain); ++c) {
        in[c][0] = in[c][1] = 0.0;
    }
    
    for (int i = 0; i < npoints; ++i) {
        int x = (int) (points[2*i  ] * dftdomain);
        int y = (int) (points[2*i+1] * dftdomain);
        fftw_complex &c = in[x + y*dftdomain];
        c[0] = 1.0;
        c[1] = 0.0;
    }
    
    fftw_plan plan;
    plan = fftw_plan_dft_2d(dftdomain, dftdomain, in, out, FFTW_FORWARD,
                            FFTW_ESTIMATE);
    fftw_execute(plan);
    
    const int offset = (dftdomain - domain) / 2;
    const int shift  = dftdomain / 2;
    
    for (int y = 0; y < domain; ++y) {
        for (int x = 0; x < domain; ++x) {
            int cx = x + offset + shift;
            int cy = y + offset + shift;
            cx = (cx >= dftdomain ? cx - dftdomain : cx);
            cy = (cy >= dftdomain ? cy - dftdomain : cy);
            fftw_complex &c = out[cx + cy*dftdomain];
            ft[2*(x + y*domain)  ] = c[0];
            ft[2*(x + y*domain)+1] = c[1];
        }
    }
    
    fftw_destroy_plan(plan);
    fftw_free(out);
    fftw_free(in);
}


void AccumPeriodogram(float *&periodogram, const float * const ft, int domain)
{
    for (int y = 0; y < domain; ++y) {
        for (int x = 0; x < domain; ++x) {
            const float &u = ft[2*(x + y*domain)  ];
            const float &v = ft[2*(x + y*domain)+1];
            periodogram[x + y*domain] += sqr(u) + sqr(v);
        }
    }
}


void RenderSpectrum(float *&image, const float * const spectrum, int domain,
                    bool mean)
{
    const float scale = 0.25f;
    
    for (int y = 0; y < domain; ++y) {
        for (int x = 0; x < domain; ++x) {
            float absolute = spectrum[x + y*domain];
            if (!mean) absolute = sqrtf(absolute);
            absolute = log2f(1.0f + scale * absolute);  // Tone mapping
            clamp1(absolute);
            image[x + y*domain] = absolute;
        }
    }
    
    // Remove DC coefficient
    image[(domain / 2) + (domain / 2) * domain] = 0.0f;
}


void SpectralStats(float *&power, float *&anisotropy, int nannuli,
                   const float * const periodogram, int domain, int npoints)
{
    assert(npoints > 1);
    
    const int domain2 = domain / 2;
    const float step  = 1.0f / nannuli;
    float *var = new float[nannuli];
    
    for (int i = 0; i < nannuli; ++i)
    {
        int rmin = (int) sqr(SQRT2 * domain2 *  i    * step);
        int rmax = (int) sqr(SQRT2 * domain2 * (i+1) * step);
        
        int Nr = 0;
        power[i] = var[i] = anisotropy[i] = 0.0f;
        
        for (int y = 0; y < domain; ++y) {
            for (int x = 0; x < domain; ++x) {
                int r = sqr(x - domain2) + sqr(y - domain2);
                if ( (r >= rmin) && (r < rmax) ) {
                    power[i] += periodogram[x + y*domain];
                    ++Nr;
                }
            }            
        }
        
        if (Nr > 0) power[i] /= (float) Nr;
        
        for (int y = 0; y < domain; ++y) {
            for (int x = 0; x < domain; ++x) {
                int r = sqr(x - domain2) + sqr(y - domain2);
                if ( (r >= rmin) && (r < rmax) ) {
                    var[i] += sqr(periodogram[x + y*domain] - power[i]);
                }
            }
        }
        
        if (Nr > 1) var[i] /= (float) (Nr - 1);
        anisotropy[i] = (sqr(power[i]) > 0.0) ? var[i] / sqr(power[i]) : 1.0f;
    }
    
    delete []var;
}

