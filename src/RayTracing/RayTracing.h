#ifndef RAY_TRACING_H
#define RAY_TRACING_H

#include <time.h>
#include <vector>
#include <algorithm>
#include <thread>
#include "../Graphics/RGB.h"
#include "../Graphics/GraphicsIO.h"
#include "../../../../Math/src/Math/Geometry/GeometricalOptics.h"
#include "Object.h"

using namespace GeometricalOptics;
using namespace ObjectLib;

#define PI 3.141592653589
#define EPS 10e-4
#define RAND_DBL (rand() / double(RAND_MAX))

/*#############################################################################
*
*						光线追踪  Ray Tracing
*
##############################################################################*/

namespace RayTracing {
	extern std::vector<Mat<>>  PointLight;													//点光源集(QuickReflect专用)
	
	extern bool haze;
	extern Mat<> haze_A;
	extern double haze_beta;
	
	extern int maxRayLevel;

	void traceRay(
		Mat<>& center, Mat<>& direct, double preSize,
		ObjectTree& objTree,
		Mat<>& R, Mat<>& G, Mat<>& B,
		int sampleSt = 0, int sampleEd = 0x7FFFFFFF
	);

	void traceRay_(
		Mat<>& center, Mat<>& direct, double width, double height,
		ObjectTree& objTree,
		Mat<>& R, Mat<>& G, Mat<>& B,
		int sampleSt = 0, int sampleEd = 0x7FFFFFFF
	);

	void traceRay_func(
		Mat<>* ScreenXVec, Mat<>* ScreenYVec, Mat<>* center, Mat<>* direct, 
		ObjectTree* objTree, double rate, 
		Mat<>* R, Mat<>* G, Mat<>* B, 
		int st, int ed
	);

	Mat<>& traceRay(ObjectTree& objTree, Mat<>& RaySt, Mat<>& Ray, Mat<>& color, int level);				//追踪光线
}

#endif