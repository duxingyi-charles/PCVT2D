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

#ifndef SPECTRAL_ANALYSIS_H
#define SPECTRAL_ANALYSIS_H


/**
 *  Computes the Fourier transform of the specified point set by inter-
 *  preting it as a superposition of Dirac peaks.
 *
 *  This is a continuous Fourier transform, accelerated via OpenMP if
 *  available.
 */
void CFT(const float * const points, int npoints, float *&ft, int domain);


/**
 *  Computes the FFT of the specified point set via FFTW.
 *
 *  The actual domain size is a multiple of the specified domain size to 
 *  increase precision. A window of the specified domain size is cut out
 *  after transformation to yield the final result.
 */
void FFT(const float * const points, int npoints, float *&ft, int domain);


/**
 *  Adds the periodogram of the specified Fourier transform to
 *  the specified periodogram.
 *
 *  The periodogram is the squared magnitude of the Fourier transform.
 */
void AccumPeriodogram(float *&periodogram, const float * const ft,
                      int domain);


/**
 *  Renders the specified spectrum to the specified float image.
 *
 *  The values of the specified spectrum are expected to be those of a
 *  mean periodogram. Thus, if 'mean' equals 'true', the values are sent
 *  unaltered to tone mapping. If, however, 'mean' equals 'false', the
 *  square root is taken before tone mapping, yielding the corresponding
 *  Fourier amplitude.
 *
 *  The resulting values are tone mapped logarithmically. The DC coefficient
 *  is removed from the image.
 */
void RenderSpectrum(float *&image, const float * const spectrum, int domain,
                    bool mean = true);


/**
 *  Computes radial power and anisotropy of the given periodogram.
 *
 *  For a good introduction to these spectral characteristics, cf.
 *  Ulichney, R., "Digital Halftoning", MIT Press, 1987.
 *
 *  For a summary related to Poisson disk distributions, cf. Lagae A. and
 *  Dutr\'{e}, P., "A Comparison of Methods for Generating Poisson Disk
 *  Distributions", Computer Graphics Forum, Vol. 27, pp. 114--129, 2008.
 */
void SpectralStats(float *&power, float *&anisotropy, int nannuli,
                   const float * const periodogram, int domain, int npoints);


#endif    // SPECTRAL_ANALYSIS_H

