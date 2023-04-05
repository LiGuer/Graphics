#ifndef GRAPHICS_ND_H
#define GRAPHICS_ND_H

#include "Graphics2D.h"
#include "Transform.h"
#include "../../../Math/src/Matrix/Operate.h"

namespace Graphics {

	/*-------------------------------- DRAW --------------------------------*/
	// Any Dim
	bool drawPoint      (Mat<ARGB>& image, Mat<int>* Z_buf, vector<int>& p);
	void drawPoint		(Mat<ARGB>& image, vector<double>& p0); 
	void drawPoint		(Mat<ARGB>& image, double x = 0, double y = 0, double z = 0); 
	void drawLine		(Mat<ARGB>& image, Mat<int>* Z_buf, vector<double>& st, vector<double>& ed);
	void drawLine		(Mat<ARGB>& image, vector<double>& st, vector<double>& ed);
	void drawSuperCuboid(Mat<ARGB>& image, vector<double>& pMin, vector<double>& pMax); 
	void drawGrid		(Mat<ARGB>& image, vector<double>& delta, vector<double>& max, vector<double>& min); 
}

#endif