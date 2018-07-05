//
// pd.cpp
//
// 2D Poisson-Disk sampling pattern generator
//
// Copyright (C) 1999, Matt Pharr <mmp@graphics.stanford.edu>
// 
// This software is placed in the public domain and is provided as is
// without express or implied warranty.
//

// Usage: pd [low] [high]
//   Generates a distribution of at least [low] but no more than [high]
//   points distributed in a posson disk distribution (i.e. no two points
//   are closer than some fixed distance together.)  A toroidal topology is
//   used, so that the pattern can be repeated across an image plane and
//   still have the Poisson Disk property at the edges.  [low] and [high]
//   should have some reasonable amount of separation  (i.e. at least 20 
//   or so)if this program is to terminate in a finite amount of time.
//
//   See Don Mitchell's classic paper, "Generating Antialiased Images at 
//   Low Sampling Densities," Proc. SIGGRAPH 1987.
//

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <cmath>
#include <cassert>
#include <ctime>

namespace DonMitchell {

#define GRID_RES 250

struct Point2 {
    double u, v;
};

// Replace these as needed with a decent random number generator.

#define INIT_RANDOM_NUMBER_GENERATOR()  srand( (unsigned)time( NULL ) );
#define UNIFORM_RANDOM()  double(rand())/double(RAND_MAX) //                drand48()

#define Point2Grid(p, r2, low, high) \
     ((low) = (int)floor(GRID_RES * ((p)-(r2))), \
      (high) = (int)ceil(GRID_RES * ((p)+(r2))))

int genDisk(int minPoints, int maxPoints, std::vector<Point2>& points) ;

int genDisk(double radius, std::vector<Point2>& points) ;

} // end of namespace

