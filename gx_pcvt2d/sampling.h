#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif


#include <Geex/mathematics/glsl_linear.h>

///
namespace Geex {

//char* SampleMode[] = {"DartThrowing","BestCandidate","Boundary","Pure","LinearPure","Penrose","Uniform"} ;

void sample(double radius, bool isTiled, bool maximize, int minMaxThrows, 
	        int multiplier, int relax, const char* method, std::vector<vec2>& points) ;

int sample_ccvt(std::vector<vec2>& points, 
		  int     NUMBER_SITES=256, 
	      double  TORUS_SIZE        = 1000,
		  bool    CONSTANT_DENSITY  = true,
		  bool    CENTROIDAL        = true) ;

void sample_Mitchell(std::vector<vec2>& points, double radius) ;

void sample_regular(std::vector<vec2>& points, double radius) ;

void sample_grid(std::vector<vec2>& points, double radius) ;

} // end of namespace Geex		