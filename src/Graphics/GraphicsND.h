#ifndef GRAPHICS_ND_H
#define GRAPHICS_ND_H

#include "../../../LiGu_Math/src/Math/Matrix/Mat.h"
#include "../RGB.h"
#include "Graphics2D.h"
#include "Transform.h"

namespace Graphics {

	/*-------------------------------- DRAW --------------------------------*/
	// Any Dim
	void drawPoint		(Mat<ARGB>& image, Mat<>& p0);				//����
	void drawPoint		(Mat<ARGB>& image, double x = 0, double y = 0, double z = 0);	//����
	void drawSuperCuboid(Mat<ARGB>& image, Mat<>& pMin, Mat<>& pMax);				//��������
	void drawGrid		(Mat<ARGB>& image, Mat<>& delta, Mat<>& max, Mat<>& min);	//������
}

#endif