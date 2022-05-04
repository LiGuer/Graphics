#ifndef INTERSECT_H
#define INTERSECT_H

#include <float.h>
#include <algorithm>
#include "../../../LiGu_Math/src/Math/Matrix/Matrix.h"

#define PI 3.141592653589

namespace Intersect {
	
extern double eps;

/*---------------- 求交点 ----------------*/
double RayPlane		(Mat<>& raySt, Mat<>& ray, Mat<>& a, double& b);						//求交-射线与平面
double RayPlane		(Mat<>& raySt, Mat<>& ray, double& A, double& B, double& C, double& D);	//求交-射线与平面(3D)
double RayCircle	(Mat<>& raySt, Mat<>& ray, Mat<>& center, double& R, Mat<>& normal);	//求交-射线与圆
double RayTriangle	(Mat<>& raySt, Mat<>& ray, Mat<>& p1, Mat<>& p2, Mat<>& p3);			//求交-射线与三角面
double RayPolygon	(Mat<>& raySt, Mat<>& ray, Mat<>* p,  int n);							//求交-射线与多边面
double RayPlaneShape(Mat<>& raySt, Mat<>& ray, Mat<>& center, Mat<>& normal, Mat<>& one, bool(*f)(double, double));//求交-射线与平面图案

double RaySphere	(Mat<>& raySt, Mat<>& ray, Mat<>& center, double& R);			//求交-射线与球面
double RaySphere	(Mat<>& raySt, Mat<>& ray, Mat<>& center, double& R, bool(*f)(double, double));			//求交-射线与球面图案
double RayCuboid	(Mat<>& raySt, Mat<>& ray, Mat<>& p1, Mat<>& p2, Mat<>& p3);			//求交-射线与矩体
double RayCuboid	(Mat<>& raySt, Mat<>& ray, Mat<>& pmin, Mat<>& pmax);					//求交-射线与矩体 (轴对齐)

}

#endif