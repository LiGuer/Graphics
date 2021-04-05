#include "RayTracing.h"
/*--------------------------------[ ��ʼ�� ]--------------------------------*/
void RayTracing::init(int width, int height) {
	g.init(width, height);
}
/*--------------------------------[ ��Ⱦ ]--------------------------------
*	[����]:
		[1] ������Ļʸ������ĻX,Y����
		[2] ����Ļÿ�����ر���
			[3] ��������ʸ��������ʸ��������׷�����
			[4] ����׷���㷨
			[5] ���ڽ�����Ƹ�����ɫ��
-------------------------------------------------------------------------*/
void RayTracing::paint() {
	//[1]
	Mat<double> ScreenVec, ScreenXVec, ScreenYVec(3, 1);
	ScreenVec.add(gCenter, Eye.negative(ScreenVec));												//��Ļ������ָ����Ļ����
	{ double t[] = { 1,-ScreenVec[0] / ScreenVec[1],0 }; ScreenYVec.getData(t).normalization(); }	//��ĻY����ʼ����Z�ᴹֱ,��z����
	ScreenXVec.crossProduct(ScreenVec, ScreenYVec).normalization();									//��ĻX��������Ļ�ᡢ��ĻY��������
	//[2]
	double minDistance = 0, RayFaceDistance;
	Mat<double> PixYVec, PixXVec, PixVec, Ray, RaySt;
	for (int x = 0; x < g.Canvas.rows; x++) {
		for (int y = 0; y < g.Canvas.cols; y++) {
			//[3]
			PixVec.add(PixXVec.mult(x - g.Canvas.rows / 2, ScreenXVec), PixYVec.mult(y - g.Canvas.cols / 2, ScreenYVec));
			Ray.add(ScreenVec, PixVec);
			RaySt.add(gCenter, PixVec);
			//[4][5]
			unsigned int color = traceRay(RaySt, Ray, 0);
			g.setPoint(x, y, color);
		}
	}
}
/*--------------------------------[ ׷�ٹ��� ]--------------------------------
*	[����]:
		[1] ���������μ����е�ÿһ��������
			[2]	�жϹ��ߺ͸��������Ƿ��ཻ�������߹����롢�������ꡢ���߼н�
			[3] ���������߹���������������ε��������
		[4] ����ù��ߵȼ�С���趨����ֵ�ȼ�
			���������η��䷽�򣬽��������Ϊ��׼���¼���
-----------------------------------------------------------------------------*/
unsigned int RayTracing::traceRay(Mat<double>& RaySt, Mat<double>& Ray, int level) {
	double minDistance = 0, RayFaceTheta, RayFaceThetaTmp;
	Mat<double> intersection, FaceVec, FaceVecTmp;
	Triangle closestTriangle;
	//[1]
	for (int i = 0; i < TriangleSet.size(); i++) {
		//[2][3]
		double distance = seekIntersection(TriangleSet[i], RaySt, Ray, FaceVecTmp, RayFaceThetaTmp, intersection);
		if (distance > 0 && distance < minDistance) {
			distance < minDistance; closestTriangle = TriangleSet[i]; FaceVec = FaceVecTmp; RayFaceTheta = RayFaceThetaTmp;
		}
	}
	//[4]
	unsigned int color = 0;
	if (closestTriangle.material != NULL && closestTriangle.material->color != 0)
		return closestTriangle.material->color;
	if (minDistance > 0 && level < maxRayLevel) {
		Mat<double> Reflect;
		Reflect.add(Ray, Reflect.add(FaceVec,Ray.negative(Reflect)));
		color = traceRay(intersection, Reflect, level + 1);
	}
	return color;
}
/*--------------------------------[ �󽻵� ]--------------------------------
*	[����]:
		[1] ����������������ʸ��
		[2] ������߾��롢��������ʸ���нǡ����������ཻ���߹��ľ���
		[3] �����������Ľ���
		[4] �жϽ����Ƿ����������ڲ�, ���񷵻�-1
*	[�㷨]:
		
---------------------------------------------------------------------------*/
double RayTracing::seekIntersection(Triangle& triangle, Mat<double>& RaySt, Mat<double>& Ray, Mat<double>& FaceVec, double& RayFaceTheta, Mat<double>& intersection) {
	//[1]
	Mat<double> edge[2], tmp;
	edge[0].add(triangle.p[1], triangle.p[0].negative(edge[0]));
	edge[1].add(triangle.p[2], triangle.p[0].negative(edge[1]));
	FaceVec.crossProduct(edge[0], edge[1]);
	//[2]
	double PointFaceDistance = fabs(FaceVec.dot(RaySt)) / (FaceVec.dot(FaceVec));
	RayFaceTheta = FaceVec.dot(Ray) / (Ray.norm() * FaceVec.norm());
	double RayFaceDistance = asin(RayFaceTheta) * PointFaceDistance;
	//[3]
	intersection.add(RaySt, intersection.mult(RayFaceDistance, Ray.normalization()));
	//[4]
	Mat<double> tmpEdge; tmpEdge.add(intersection, triangle.p[0].negative(tmpEdge));
	double inverDeno = 1 / (edge[0].dot(edge[0]) * edge[1].dot(edge[1]) - edge[0].dot(edge[0]) * edge[0].dot(edge[0]));
	double u = (edge[1].dot(edge[1]) * edge[0].dot(tmpEdge) - edge[0].dot(edge[1]) * edge[1].dot(tmpEdge)) * inverDeno;
	if (u < 0 || u > 1)  return -1;// if u out of range, return directly
	double v = (edge[0].dot(edge[0]) * edge[1].dot(tmpEdge) - edge[0].dot(edge[1]) * edge[0].dot(tmpEdge)) * inverDeno;
	if (v < 0 || v > 1)  return -1;// if v out of range, return directly
	return u + v <= 1 ? RayFaceDistance : -1;
}