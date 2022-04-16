#ifndef GRAPHICS_TRANSFORM_H
#define GRAPHICS_TRANSFORM_H

#include "../../../LiGu_Math/src/Math/Matrix/Transform.h"
#include "../../../LiGu_Math/src/Math/Matrix/Submatrix.h"
#include "../RGB.h"
#include "Graphics2D.h"

namespace Graphics {
	extern Mat<> TransformMat;

	/*---------------- ���α任 Transformation ----------------*/
	Mat<>& translate(Mat<>& delta);						//ƽ��
	Mat<>& rotate	(Mat<>& theta, Mat<>& center);		//��ת
	Mat<>& rotate	(Mat<>& rotateAxis, double theta, Mat<>& center);	//��ת 3D
	Mat<>& scale	(Mat<>& ratio,Mat<>& center);		//����
	Mat<>& reflect	(Mat<>& e, Mat<>& center);			//����
}

#endif