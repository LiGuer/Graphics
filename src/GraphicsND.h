#ifndef GRAPHICS_ND_H
#define GRAPHICS_ND_H

#include "Graphics2D.h"
#include "Transform.h"
#include "../../../../Math/src/Matrix/Operate.h"

namespace Graphics {

	/*-------------------------------- DRAW --------------------------------*/
	// Any Dim
	bool drawPoint      (Mat<ARGB>& image, Mat<int>* Z_buf, Mat<int>& p) 
	void drawPoint		(Mat<ARGB>& image, Mat<>& p0); 
	void drawPoint		(Mat<ARGB>& image, double x = 0, double y = 0, double z = 0); 
	void drawLine		(Mat<ARGB>& image, Mat<int>* Z_buf, Mat<>& st, Mat<>& ed);
	void drawLine		(Mat<ARGB>& image, Mat<>& st, Mat<>& ed);
	void drawSuperCuboid(Mat<ARGB>& image, Mat<>& pMin, Mat<>& pMax); 
	void drawGrid		(Mat<ARGB>& image, Mat<>& delta, Mat<>& max, Mat<>& min); 
}

#endif