#include "Transform.h"

Mat<> Graphics::TransformMat(4);

/*#############################################################################

*                    Transformation 

##############################################################################*/

// ƽ��
Mat<>& Graphics::translate(Mat<>& delta) {
	static Mat<> translateMat;
	Matrix::E(translateMat.alloc(TransformMat.rows, TransformMat.cols));

	for (int i = 0; i < delta.rows; i++) 
		translateMat(i + 1, 0) = delta[i];

	Matrix::mul(TransformMat, translateMat, TransformMat);
	return TransformMat;
}

// ��ת
Mat<>& Graphics::rotate(Mat<>& theta, Mat<>& center) {
	static Mat<> tmp, rotateMat;
	translate(Matrix::negative(tmp, center));

	Matrix::rotate(theta, rotateMat);
	Matrix::E(tmp.alloc(TransformMat.rows, TransformMat.cols));
	Matrix::setBlock(tmp, rotateMat, 1, 1);

	Matrix::mul(TransformMat, tmp, TransformMat);
	translate(center);
	return TransformMat;
}

// 3D����Ԫ��
Mat<>& Graphics::rotate(Mat<>& rotateAxis, double theta, Mat<>& center) {
	if (TransformMat.rows != 4 || TransformMat.cols != 4)
		exit(-1);

	static Mat<> tmp, rotateMat;
	translate(Matrix::negative(tmp, center));

	Matrix::rotate(rotateAxis, theta, rotateMat);

	Matrix::mul(TransformMat, rotateMat, TransformMat);
	translate(center);
	return TransformMat;
}

// ����
Mat<>& Graphics::scale(Mat<>& ratio, Mat<>& center) {
	static Mat<> tmp, scaleMat;
	translate(Matrix::negative(tmp, center));

	Matrix::scale(ratio, scaleMat);
	Matrix::E(tmp.alloc(TransformMat.rows, TransformMat.cols));
	Matrix::setBlock(tmp, scaleMat, 1, 1);

	Matrix::mul(TransformMat, tmp, TransformMat);
	translate(center);
	return TransformMat;
}

// ����
Mat<>& Graphics::reflect(Mat<>& e, Mat<>& center) {
	Mat<> tmp, reflectMat;
	translate(Matrix::negative(tmp, center));

	Matrix::reflect(e, reflectMat);
	Matrix::E(tmp.alloc(TransformMat.rows, TransformMat.cols));
	Matrix::setBlock(tmp, reflectMat, 1, 1);

	Matrix::mul(TransformMat, tmp, TransformMat);
	translate(center);
	return TransformMat;
}

// ͸��ͶӰ
/*
Mat<>& Graphics::perspect(Mat<>& e, Mat<>& center) {

}*/