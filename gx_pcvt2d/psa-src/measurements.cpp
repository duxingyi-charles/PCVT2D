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

#include "measurements.h"

#include "util.h"


void BasicStats(const std::list<float> &values, struct stats &s)
{
    double m_old = 0, m_new = 0, m2 = 0;
    double delta, sd;
    
    s.min = FLT_MAX;
    s.max = FLT_MIN;
    
    std::list<float>::const_iterator val;
    int i = 0;
    for (val = values.begin(); val != values.end(); ++val) {
        delta = *val - m_old;
        m_new = m_old + delta / (i + 1.0);
        m2    = m2 + delta * (*val - m_new);
        m_old = m_new;
        s.min = std::min(s.min, *val);
        s.max = std::max(s.max, *val);
        ++i;
    }
    sd = sqrt(m2 / (double) values.size());
    
    s.mean = (float) m_new;
    s.sd   = (float) sd;
}


float L2normStarDiscrepancy(const float * const points, int npoints)
{
    double s1 = 0, s2 = 0, s3 = 0, s4 = 0;
    
    for (int i = 0; i < npoints; ++i) {
        for (int j = i + 1; (j < npoints) && (i != npoints - 1); ++j) {
            s1 += (1.0 - std::max(points[2*i  ], points[2*j]  )) *
                  (1.0 - std::max(points[2*i+1], points[2*j+1]));
        }
        s2 += (1.0 - points[2*i]) * (1.0 - points[2*i+1]);
        s3 += (1.0 - sqr(points[2*i])) * (1.0 - sqr(points[2*i+1]));
    }
    s1 *= 2.0;
    s3 *= npoints / 2.0;
    s4  = sqr(npoints) / 9.0;
    
    return (float) (sqrt(s1 + s2 - s3 + s4) / (double) npoints);
}


void DistanceMeasures(const float * const points, int npoints,
                      float &mindist, float &avgmindist)
{
    mindist = FLT_MAX;
    avgmindist = 0;
    
    for (int i = 0; i < npoints; ++i) {
        float local_md = FLT_MAX;
        for (int j = 0; j < npoints; ++j) {
            if (i == j) continue;
            float d = dist_sqr_torus(points[2*i], points[2*i+1],
                                     points[2*j], points[2*j+1]);
            local_md = std::min(d, local_md);
        }
        mindist = std::min(local_md, mindist);
        avgmindist += sqrt(local_md);
    }
    
    mindist = sqrt(mindist);
    avgmindist /= npoints;
}

