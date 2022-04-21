#ifndef RAY_TRACING_H
#define RAY_TRACING_H

#include <time.h>
#include <vector>
#include <algorithm>
#include "../RGB.h"
#include "../GraphicsIO.h"
#include "GeometricalOptics.h"
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
	/*
	bool haze = false;
	Mat<> haze_A{ 3 };
	double haze_beta = 1;*/
	extern int maxRayLevel;

	void traceRay(
		Mat<>& center, Mat<>& direct, double width, double height,
		ObjectTree& objTree,
		Mat<>& R, Mat<>& G, Mat<>& B,
		int sampleSt = 0, int sampleEd = 0x7FFFFFFF
	);

	Mat<>& traceRay(ObjectTree& objTree, Mat<>& RaySt, Mat<>& Ray, Mat<>& color, int level);				//追踪光线
}

#endif