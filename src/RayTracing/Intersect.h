#ifndef INTERSECT_H
#define INTERSECT_H

#include <float.h>
#include <algorithm>
#include "../../../LiGu_Math/src/Math/Matrix/Matrix.h"

#define PI 3.141592653589
#define EPS 10e-4

namespace Intersect {
	
/*---------------- 求交点 ----------------*/
double RayPlane		(Mat<>& RaySt, Mat<>& Ray, double& A, double& B, double& C, double& D);	//求交-射线与平面
double RayCircle	(Mat<>& RaySt, Mat<>& Ray, Mat<>& Center, double& R, Mat<>& normal);	//求交-射线与圆
double RayTriangle	(Mat<>& RaySt, Mat<>& Ray, Mat<>& p1, Mat<>& p2, Mat<>& p3);	//求交-射线与三角面
double RayPolygon	(Mat<>& RaySt, Mat<>& Ray, Mat<>* p,  int n);					//求交-射线与多边面
double RayPlaneShape(Mat<>& RaySt, Mat<>& Ray, Mat<>& Center, Mat<>& normal, Mat<>& one, bool(*f)(double, double));//求交-射线与平面图形
double RaySphere	(Mat<>& RaySt, Mat<>& Ray, Mat<>& center, double& R, bool(*f)(double, double) = NULL);			//求交-射线与球
double RayCuboid	(Mat<>& RaySt, Mat<>& Ray, Mat<>& p1, Mat<>& p2, Mat<>& p3);	//求交-射线与长方体
double RayCuboid	(Mat<>& RaySt, Mat<>& Ray, Mat<>& pmin, Mat<>& pmax);			//求交-射线与长方体 (轴对齐)

}

#endif