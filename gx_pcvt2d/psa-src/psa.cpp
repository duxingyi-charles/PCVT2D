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
 
#include <iostream>
#include <cstdio>
#include <list>

#include "fileio.h"
#include "measurements.h"
#include "spectralanalysis.h"
#include "util.h"
#include "psa.h"

using namespace std;

namespace Geex {

void usage()
{
    cout << "usage: psa [-fes] <file or folder>\n"
            "  [h]: show this message\n"
            "  [f]: use a FFT instead of a continuous FT\n"
            "  [e]: output EPS instead of PDF files\n"
            "  [s]: output summary instead of separate files\n";
}

//int main (int argc, char * const argv[])
//{
//    const string pname = "psa";
//    string fname;
//    string ext;
//    
//    bool fft = false;
//    bool eps = false;
//    bool summary = false;
//    
//    
//    // Query arguments
//    if (argc < 2 || argc > 3) {
//        usage(); exit(0);
//    } 
//    else
//    {
//        const bool hasOption = (argv[1][0] == '-');
//        
//        if (hasOption) {
//            string options = string(argv[1]);
//            for (unsigned i = 1; i < options.length(); ++i) {
//                switch (options[i]) {
//                    case 'h': usage(); exit(0);
//                    case 'f': fft = true; break;
//                    case 'e': eps = true; break;
//                    case 's': summary = true; break;
//                    default:
//                        cout << pname << ": illegal option " << options[i]
//                             << endl;
//                        usage(); exit(0);
//                }
//            }
//            if (argc != 3) {
//                cout << pname << ": no file or directory given" << endl;
//                usage(); exit(0);
//            }
//        }
//        fname = string(argv[argc-1]);
//        ext   = eps ? ".eps" : ".pdf";
//    }
    
int psa (bool fft, bool eps, bool summary, std::string& fname) {    
    const string pname = "psa";
    string ext = eps ? ".eps" : ".pdf";

    // Setup a list of all file names to be processed. This list consists of
    // either the single specified file or all files in the specified folder.
    list<string> fnames;
    FileList(fname, fnames);
    
    if (fnames.size() == 0) {
        if ( FolderExists(fname) )
            cerr << "No " << FILE_EXT_RPS << " or " << FILE_EXT_TXT
            << " files found in " << fname << endl;
        else
            cerr << "No " << FILE_EXT_RPS << " or " << FILE_EXT_TXT
            << " file " << fname << endl;
        //exit(1);
		return 1 ;
    }
    
    
    // Read the number of points from the first file to determine boundary
    // conditions, e.g. fourier domain.
    int npoints = 0;
    if ( !ReadNumPoints(fnames.front(), npoints) ) {
        //exit(1);
		return 1; 
    }

    // Perform the analysis
    float *points = 0;
    int nsets = 0;
    int lastnpoints = npoints;
    
    const double d_max   = sqrt(2.0 / (SQRT3 * npoints));
    const int ftdomain   = nextpow2(2 * (int)(4.0 / d_max));
    const int nannuli    = SQRT2 * ftdomain / 2;
    
    float *ft            = new float[ftdomain * ftdomain * 2];
    float *periodogram   = new float[ftdomain * ftdomain];
    float *spectrumImage = new float[ftdomain * ftdomain];
    float *power         = new float[nannuli];
    float *anisotropy    = new float[nannuli];
    
    SetZero(ft, ftdomain * ftdomain * 2);
    SetZero(periodogram, ftdomain * ftdomain);
    
    list<float> l2disc, mindist, mindist_norm, avgmd_norm;
    stats s_l2disc, s_mindist, s_mindist_norm, s_avgmd_norm;

    
    list<string>::iterator file;
    for (file = fnames.begin(); file != fnames.end(); ++file)
    {
        cout << "Processing " << *file << endl;
        
        if (points) delete []points;
        
        if (!LoadPoints(*file, &points, npoints))
            continue;
        
        if (npoints != lastnpoints) {
            cerr << "Processed files have differing number of points" << endl;
            cerr << *file << " has " << npoints << " points while previous "
                 << "had " << lastnpoints << endl;
            //exit(1);
			//return 1;
        }
        
        // Fourier transform
        if (fft)
            FFT(points, npoints, ft, ftdomain);
        else
            CFT(points, npoints, ft, ftdomain);
        
        AccumPeriodogram(periodogram, ft, ftdomain);
        
        // Measurements
        float md, amd;
        DistanceMeasures( points, npoints, md, amd );
        l2disc.push_back( L2normStarDiscrepancy(points, npoints) );
        mindist.push_back(md);
        mindist_norm.push_back(md / d_max);
        avgmd_norm.push_back(amd / d_max);
        
        ++nsets;
        lastnpoints = npoints;
    }
    
    if (nsets == 0) {
        cerr << "No valid point sets" << endl;
        //exit(1);
		return 1 ;
    }
    
    // Normalization
    Divide(periodogram, ftdomain * ftdomain, npoints * nsets);
    
    
    // Compute radial power and anisotropy
    SpectralStats(power, anisotropy, nannuli, periodogram, ftdomain, npoints);
    Decibel(anisotropy, nannuli);

    // Render spectrum. In case of multiple files, we render the mean
    // periodogram; for single point sets, we render the amplitude spectrum.
    RenderSpectrum(spectrumImage, periodogram, ftdomain, (nsets > 1));
    
    // Basic statistical analysis for the measurements:
    // mean, std. dev., min and max value
    BasicStats(l2disc, s_l2disc);
    BasicStats(mindist, s_mindist);
    BasicStats(mindist_norm, s_mindist_norm);
    BasicStats(avgmd_norm, s_avgmd_norm);
    
    
    // Generate output
    if (FolderExists(fname)) {
        RemoveTrailingSlash(fname);
    }
    string base = GetFileName(fname);
    
    if (summary)
    {
        string c_fname = base + ext;              // Single composite
        SaveComposite(c_fname, eps, points, npoints, spectrumImage, ftdomain,
                      ftdomain, power, anisotropy, nannuli, s_l2disc,
                      s_mindist, s_mindist_norm, s_avgmd_norm,
                      nsets, 512, 512);
    }
    else
    {
        string p_fname = base + "_pts" + ext;     // (Last) point set
        SavePointPlot(p_fname, eps, points, npoints, 512, 512);
        
        string png_fname = base + "_spec.png";    // Spectrum
        SaveFloatPNG(png_fname, spectrumImage, ftdomain, ftdomain);

        string rp_fname = base + "_rp" + ext;     // Radial power
        SaveSpectralPlot(rp_fname, eps, "power", power, nannuli, 512, 320,
                         info_rp);
        
        string ani_fname = base + "_ani" + ext;   // Anisotropy
        SaveSpectralPlot(ani_fname, eps, "anisotropy", anisotropy, nannuli,
                         512, 320, ani_format(nsets));
    }
    
    // Print measurements
    if (nsets == 1) {
        printf("L2-norm discr: %.6f\n", s_l2disc.mean);
        printf("Mindist      : %.6f\n", s_mindist.mean);
        printf("Norm. mindist: %.6f\n", s_mindist_norm.mean);
        printf("Norm. avg. md: %.6f\n", s_avgmd_norm.mean);
    } else {
        printf("L2-norm discr: mean=%.6f sd=%.6f min=%.6f max=%.6f\n",
               s_l2disc.mean, s_l2disc.sd, s_l2disc.min, s_l2disc.max);
        printf("Mindist      : mean=%.6f sd=%.6f min=%.6f max=%.6f\n",
               s_mindist.mean, s_mindist.sd, s_mindist.min, s_mindist.max);
        printf("Norm. mindist: mean=%.6f sd=%.6f min=%.6f max=%.6f\n",
               s_mindist_norm.mean, s_mindist_norm.sd, s_mindist_norm.min,
               s_mindist_norm.max);
        printf("Norm. avg. md: mean=%.6f sd=%.6f min=%.6f max=%.6f\n",
               s_avgmd_norm.mean, s_avgmd_norm.sd, s_avgmd_norm.min,
               s_avgmd_norm.max);
    }

    delete []ft;
    delete []periodogram;
    delete []spectrumImage;
    delete []power;
    delete []anisotropy;
    return 0;
}

}
