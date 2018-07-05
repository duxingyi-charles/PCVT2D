#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#include <algorithm>
#include <cmath>
#include <cstring>
#include <ctime>
#include <fstream>

#include "sampling.h"
#include "pd.h"
#include "PDSample/PDSampling.h"
#include "ccvt/ccvt_metric.h"
#include "ccvt/ccvt_optimizer.h"
#include "ccvt/ccvt_point.h"
#include "ccvt/ccvt_site.h"


///
namespace Geex {

//void sample_regular(std::vector<vec2>& points, double radius) {
//}
//
//void sample_grid(std::vector<vec2>& points, double radius) {
//	int    n = ceil(1.0/radius) ;
//	double len = 1.0/n ;
//	for(int i=0; i<n; ++i) {
//		for(int j=0; j<n; ++j) {
//			points.push_back(vec2((i+0.5)*len, (j+0.5)*len)) ;
//		}
//	}
//}
//
//void sample_Mitchell(std::vector<vec2>& points, double radius) {
//	std::vector<DonMitchell::Point2> disks ;
//	DonMitchell::genDisk(radius, disks) ;
//	for(unsigned int i=0; i<disks.size(); ++i) {
//		points.push_back(vec2(disks[i].u, disks[i].v)) ;
//	}
//}
//
//// discrete space with constant density;
//// the points form a regular grid
//void constant_regular_density(ccvt::Point2::List& points, const int numberOfPoints, const double torusSize) {
//  double n = sqrt(static_cast<double>(numberOfPoints));
//  for (int x = 0; x < n; ++x) {
//    for (int y = 0; y < n; ++y) {
//      double dx = x / n * torusSize;
//      double dy = y / n * torusSize;
//      points.push_back(ccvt::Point2(dx, dy));
//    }
//  }
//}
//
//// discrete space with constant density;
//// the points are randomly distributed
//void constant_random_density(ccvt::Point2::List& points, const int numberOfPoints, const double torusSize) {
//  for (int i = 0; i < numberOfPoints; ++i) {
//    double x = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * torusSize;
//    double y = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * torusSize;
//    points.push_back(ccvt::Point2(x, y));
//  }
//}
//
//// discrete space with the density function e^(-20x^2-20y^2)+0.2sin^2(PIx)sin^2(PIy);
//// the points are generated via rejection sampling
//void nonconstant_density(ccvt::Point2::List& points, const int numberOfPoints, const double torusSize) {
//  const double E = 2.718281828459;
//  const double PI = 3.141592653590;
//  while (points.size() < static_cast<unsigned int>(numberOfPoints)) {
//    double x = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * 2 - 1;
//    double y = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * 2 - 1;
//    double p = pow(E, -20.0 * x * x - 20.0 * y * y) + 0.2 * sin(PI * x) * sin(PI * x) * sin(PI * y) * sin(PI * y);
//    double r = static_cast<double>(rand() % RAND_MAX) / RAND_MAX;
//    if (p >= r) {
//      points.push_back(ccvt::Point2((x + 1) / 2 * torusSize, (y + 1) / 2 * torusSize));
//    }
//  }
//}
//
//// export sites to an EPS image
//bool save_eps(const char* filename, const ccvt::Site<ccvt::Point2>::Vector& sites, const double width, const double height, const double radius) {
//  std::ofstream stream(filename, std::ios::out);
//  if (stream.bad()) {
//    return false;
//  }
//
//  stream << "%!PS-Adobe EPSF-3.0\n";
//  stream << "%%HiResBoundingBox: " << 0.0 << " " << 0.0 << " " << width << " " << height << "\n";
//  stream << "%%BoundingBox: " << 0 << " " << 0 << " " << static_cast<int>(width) << " " << static_cast<int>(height) << "\n";
//  stream << "\n";
//  stream << "%% Sites: " << sites.size() << "\n";
//  stream << "\n";
//  stream << "/radius { " << radius << " } def\n";
//  stream << "\n";
//  stream << "/p { radius 0 360 arc closepath fill } def\n";
//  stream << "\n";
//  stream << "0 0 0 setrgbcolor\n";
//  stream << "\n";
//  for (unsigned int i = 0; i < sites.size(); ++i) {
//    stream << sites[i].location.x << " " << sites[i].location.y << " p\n";
//  }
//  stream << "\n";
//  stream << "showpage\n";
//
//  stream.close();
//  return true;
//}
//
//int sample_ccvt(std::vector<vec2>& samples, int     NUMBER_SITES, 
//	      double  TORUS_SIZE,
//		  bool    CONSTANT_DENSITY,
//		  bool    CENTROIDAL) {
//typedef ccvt::Optimizer<ccvt::Site<ccvt::Point2>, ccvt::Point2, ccvt::MetricToroidalEuclidean2> Optimizer;
//
//  int NUMBER_POINTS = 1024*NUMBER_SITES;
//  // intializing the underlying discrete space
//  ccvt::Point2::List points;
//  if (CONSTANT_DENSITY) {
//    constant_regular_density(points, NUMBER_POINTS, TORUS_SIZE);
//  } else {
//    nonconstant_density(points, NUMBER_POINTS, TORUS_SIZE);
//  }
//
//  // initializing the Voronoi sites with equal capacity
//  unsigned int overallCapacity = static_cast<int>(points.size());
//  ccvt::Site<ccvt::Point2>::List sites;
//  for (int i = 0; i < NUMBER_SITES; ++i) {
//    double x = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * TORUS_SIZE;
//    double y = static_cast<double>(rand() % RAND_MAX) / RAND_MAX * TORUS_SIZE;
//    int capacity = overallCapacity / (NUMBER_SITES - i);
//    overallCapacity -= capacity;
//    sites.push_back(ccvt::Site<ccvt::Point2>(i, capacity, ccvt::Point2(x, y)));
//  }
//
//  clock_t start = clock();
//
//  // initializing the CCVT
//  clock_t startInitialization = clock();
//  printf("initialization...");
//  Optimizer optimizer;
//  ccvt::MetricToroidalEuclidean2 metric(ccvt::Point2(TORUS_SIZE, TORUS_SIZE));
//  optimizer.initialize(sites, points, metric);
//  printf("done\n");
//  clock_t endInitialization = clock();
//
//  // optimization
//  int iteration = 0;
//  bool stable;
//  do {
//    printf("iteration %d...", ++iteration);
//    stable = optimizer.optimize(CENTROIDAL);
//    printf("done\n");
//  } while (!stable);
//  
//  clock_t end = clock();
//
//  const ccvt::Site<ccvt::Point2>::Vector& result = optimizer.sites();
//  
//  //// writing the Voronoi sites to console
//  //if (RESULT_PRINT) {
//  //  printf("\nresult:\n");
//  //  for (unsigned int i = 0; i < result.size(); ++i) {
//  //    printf("site %d: %f, %f\n", result[i].id, result[i].location.x, result[i].location.y);
//  //  }
//  //}
//
//  printf("\ninitialization time: %.3f sec\n", static_cast<double>(endInitialization - startInitialization) / CLOCKS_PER_SEC);
//  printf("computation time: %.3f sec\n", static_cast<double>(end - start) / CLOCKS_PER_SEC);
//
//  // writing the Voronoi sites to EPS file
//  //if (RESULT_FILE) {
//  //  if (save_eps(RESULT_FILENAME, result, TORUS_SIZE, TORUS_SIZE, RESULT_RADIUS)) {
//  //    printf("\nresult saved in '%s'\n", RESULT_FILENAME);
//  //  } else {
//  //    printf("\nresult could not be saved in '%s'\n", RESULT_FILENAME);
//  //  }
//  //}
//
//    for (unsigned int i = 0; i < result.size(); ++i) {
//		samples.push_back(vec2((double)(result[i].location.x)/TORUS_SIZE, (double)(result[i].location.y)/TORUS_SIZE));
//    }
//
//  printf("\n");
//
//  return 0;	
//}
//
//
//void usage(char *app)
//{
//	printf(	"Usage: %s [-m] [-t] [-r <relax count=0>] [-M <multiplier=10>] [-N <minMaxThrows=1000>] <method> <radius> <output>\n", app);
//	printf(	"	-m		maximize point set after sampling\n");
//	printf(	"	-t		use tiled domain\n");
//	printf(	"	-r		apply the specified number of relaxations after sampling (requires qvoronoi)\n");
//	printf(	"	-M		set multiplier for DartThrowing and BestCandidate samplers\n");
//	printf(	"	-N		set minimum number of maximum throws for DartThrowing sampler\n");
//	printf(	"	available methods = {\n"
//			"		DartThrowing, (uses multiplier and minMaxThrows)\n"
//			"		BestCandidate, (uses multiplier)\n"
//			"		Boundary, \n"
//			"		Pure, \n"
//			"		LinearPure, \n"
//			"		Penrose, \n"
//			"		Uniform" "}\n");
//
//	exit(1);
//}

void sample(double radius, bool isTiled, bool maximize, int minMaxThrows, 
	        int multiplier, int relax, const char* method, std::vector<vec2>& points) {
	//PDSampler *sampler;
	//FILE *output;
	//double startTime;
	//float elapTime;
	//int i, N;

	//if (radius<0.0005 || radius>.2) {
	//	printf("Radius (%f) is outside allowable range.\n", radius);
	//	exit(1);
	//}

	//if (!strcmp(method, "DartThrowing")) {
	//	sampler = new DartThrowing(radius, isTiled, minMaxThrows, multiplier);
	//} else if (!strcmp(method, "BestCandidate")) {
	//	sampler = new BestCandidate(radius, isTiled, multiplier);
	//} else if (!strcmp(method, "Boundary")) {
	//	sampler = new BoundarySampler(radius, isTiled);
	//} else if (!strcmp(method, "Pure")) {
	//	if (!isTiled) {
	//		printf("Pure sampler does not support untiled domain.\n");
	//		return ; //exit(1);
	//	}
	//	sampler = new PureSampler(radius);
	//} else if (!strcmp(method, "LinearPure")) {
	//	if (!isTiled) {
	//		printf("LinearPure sampler does not support untiled domain.\n");
	//		return ; //exit(1);
	//	}
	//	sampler = new LinearPureSampler(radius);
	//} else if (!strcmp(method, "Penrose")) {
	//	if (isTiled) {
	//		printf("Penrose sampler does not support tiled domain.\n");
	//		return ; //exit(1);
	//	}
	//	sampler = new PenroseSampler(radius);
	//} else if (!strcmp(method, "Uniform")) {
	//	sampler = new UniformSampler(radius);
	//} else {
	//	//printf("Unrecognized sampling method (%s).\n", );
	//	return ; //exit(1);
	//}

	//startTime = timeInSeconds();
	//sampler->complete();
	//if (maximize) sampler->maximize();
	//for (i=0; i<relax; i++) sampler->relax();
	//elapTime = (float) (timeInSeconds() - startTime);

	//N = (int) sampler->points.size();
	//for (i=0; i<N; i++) {
	//	Vec2 pt = sampler->points[i] ;
	//	points.push_back(vec2(pt.x, pt.y));
	//}

	////output = fopen(outputPath,"wb");
	////if (!output) {
	////	printf("Unable to open output file (%s).\n", outputPath);
	////	exit(1);
	////} 
	//std::cout << "number of samples: " << N << std::endl ;
	//std::cout << "elapped time: " << elapTime << std::endl ;
	//std::cout << "sample radius: " << radius << std::endl ;

	////fwrite(&N, 4, 1, output);
	////fwrite(&elapTime, 4, 1, output);
	////fwrite(&radius, 4, 1, output);
	////fwrite(&isTiled, 4, 1, output);
	////for (i=0; i<N; i++) {
	////	fwrite(&sampler->points[i], 8, 1, output);
	////}
	////fclose(output);
	////printf("Wrote: %s (%d points, %fs, %6.1f pts/s)\n", outputPath, N, elapTime, N/elapTime);

	//delete sampler;
}

} // end of namespace Geex