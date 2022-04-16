#ifndef INTERSECT_H
#define INTERSECT_H

#include "../../../LiGu_Math/src/Math/Matrix/Matrix.h"

namespace Intersect {
	
/*---------------- �󽻵� ----------------*/
double RayPlane		(Mat<>& RaySt, Mat<>& Ray, double& A, double& B, double& C, double& D);	//��-������ƽ��
double RayCircle	(Mat<>& RaySt, Mat<>& Ray, Mat<>& Center, double& R, Mat<>& normal);	//��-������Բ
double RayTriangle	(Mat<>& RaySt, Mat<>& Ray, Mat<>& p1, Mat<>& p2, Mat<>& p3);	//��-������������
double RayPolygon	(Mat<>& RaySt, Mat<>& Ray, Mat<>* p,  int n);					//��-����������
double RayPlaneShape(Mat<>& RaySt, Mat<>& Ray, Mat<>& Center, Mat<>& normal, Mat<>& one, bool(*f)(double, double));//��-������ƽ��ͼ��
double RaySphere	(Mat<>& RaySt, Mat<>& Ray, Mat<>& center, double& R, bool(*f)(double, double) = NULL);			//��-��������
double RayCuboid	(Mat<>& RaySt, Mat<>& Ray, Mat<>& p1, Mat<>& p2, Mat<>& p3);	//��-�����볤����
double RayCuboid	(Mat<>& RaySt, Mat<>& Ray, Mat<>& pmin, Mat<>& pmax);			//��-�����볤���� (�����)

}

#endif