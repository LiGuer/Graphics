#ifndef RAY_TRACING_MATERIAL_H
#define RAY_TRACING_MATERIAL_H

#include "../../../LiGu_Math/src/Math/Matrix/Matrix.h"

struct Material {															//²ÄÖÊ
	Mat<> color{ 3 }, refractRate{ 3 };
	bool
		rediate = 0,
		quickReflect = 0,
		diffuseReflect = 0;
	double
		reflect = 0, reflectLossRate = 1,
		refract = 0, refractLossRate = 1;
};

#endif