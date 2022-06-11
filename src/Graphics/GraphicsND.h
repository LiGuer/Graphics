#ifndef GRAPHICS_ND_H
#define GRAPHICS_ND_H

#include "../../../LiGu_Math/src/Math/Matrix/Mat.h"
#include "../RGB.h"
#include "Graphics2D.h"
#include "Transform.h"

namespace Graphics {

	/*-------------------------------- DRAW --------------------------------*/
	// Any Dim
	void drawPoint		(Mat<ARGB>& image, Mat<>& p0);				//画点
	void drawPoint		(Mat<ARGB>& image, double x = 0, double y = 0, double z = 0);	//画点
	void drawLine		(Mat<ARGB>& image, Mat<>& st, Mat<>& ed);
	void drawSuperCuboid(Mat<ARGB>& image, Mat<>& pMin, Mat<>& pMax);				//画立方体
	void drawGrid		(Mat<ARGB>& image, Mat<>& delta, Mat<>& max, Mat<>& min);	//画网格
}

#endif