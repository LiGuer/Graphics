#ifndef GEOMETERICAL_OPTICS_H
#define GEOMETERICAL_OPTICS_H

#include "../../../LiGu_Math/src/Math/Matrix/Matrix.h"
using namespace Matrix;

#define PI 3.141592653589
#define RAND_DBL (rand() / double(RAND_MAX))

namespace GeometricalOptics {

	/*---- ���� ----*/
	inline Mat<>& reflect(Mat<>& RayI, Mat<>& faceVec, Mat<>& RayO) {
		mul(RayO, -2 * dot(faceVec, RayI), faceVec);
		add(RayO, RayO, RayI);
		return normalize(RayO);
	}

	/*---- ���� ----*/
	inline Mat<>& refract(Mat<>& RayI, Mat<>& faceVec, Mat<>& RayO, double rateI, double rateO) {
		double k = rateI / rateO,
			CosI = dot(faceVec, RayI),
			CosO = 1 - pow(k, 2) * (1 - pow(CosI, 2));

		if (CosO < 0)				//ȫ����
			return reflect(RayI, faceVec, RayO);

		mul(RayO, -CosI - (CosI > 0 ? -1 : 1) * sqrt(CosO) / k, faceVec);
		add(RayO, RayO, RayI);
		return normalize(RayO);
	}

	/*---- ������ ----*/
	inline Mat<>& diffuseReflect(Mat<>& RayI, Mat<>& faceVec, Mat<>& RayO) {
		double r1 = 2 * PI * RAND_DBL, r2 = RAND_DBL;
		static Mat<> t(3), u, v;
		mul(faceVec, dot(faceVec, RayI) > 0 ? -1 : 1, faceVec);
		t[0] = fabs(faceVec[0]) > 0.1 ? 0 : 1;
		t[1] = t[0] == 0 ? 1 : 0;
		mul(u, cos(r1) * sqrt(r2), normalize(cross_(u, t, faceVec)));
		mul(v, sin(r1) * sqrt(r2), normalize(cross_(v, faceVec, u)));

		mul(RayO, sqrt(1 - r2), faceVec);
		add(RayO, RayO, mul(u, u, v));
		return normalize(RayO);
	}

	/*---- �� (����ͬ��) ----*/
	inline double Haze(double I, double A, double dis, double beta) {
		double t = exp(-beta * dis);
		return I * t + A * (1 - t);
	}

	inline Mat<>& Haze(Mat<>& I, Mat<>& O, Mat<>& A, double dis, double beta) {
		static Mat<> tmp(3);
		double t = exp(-beta * dis);
		mul(O, t, I);
		mul(O, O, mul(tmp, 1 - t, A));
		return O;
	}
}

#endif