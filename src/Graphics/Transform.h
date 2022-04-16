#ifndef GRAPHICS_TRANSFORM_H
#define GRAPHICS_TRANSFORM_H

#include "../../../LiGu_Math/src/Math/Matrix/Transform.h"
#include "../../../LiGu_Math/src/Math/Matrix/Submatrix.h"
#include "../RGB.h"
#include "Graphics2D.h"

namespace Graphics {
	extern Mat<> TransformMat;

	/*---------------- 几何变换 Transformation ----------------*/
	Mat<>& translate(Mat<>& delta);						//平移
	Mat<>& rotate	(Mat<>& theta, Mat<>& center);		//旋转
	Mat<>& rotate	(Mat<>& rotateAxis, double theta, Mat<>& center);	//旋转 3D
	Mat<>& scale	(Mat<>& ratio,Mat<>& center);		//缩放
	Mat<>& reflect	(Mat<>& e, Mat<>& center);			//镜像
}

#endif