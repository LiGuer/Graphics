#include "RayTracing.h"
/*--------------------------------[ ��ʼ�� ]--------------------------------*/
void RayTracing::init(int width, int height) {
	g.init(width, height);
}
/*--------------------------------[ ��Ⱦ ]--------------------------------*/
void RayTracing::paint() {
	double minDistance = 0, RayFaceDistance;
	
	Mat<double> ScreenVec, PixYVec(3, 1), PixXVec, Ray, RaySt;
	ScreenVec.add(gCenter, ScreenVec.negative(Eye));
	for (int x = 0; x < g.Canvas.rows; x++) {
		for (int y = 0; y < g.Canvas.cols; y++) {
			//Compute Ray
			double t[] = { 1,-ScreenVec[0] / ScreenVec[1],0 };
			PixYVec.mult(y, PixYVec.getData(t));								//��ĻY����ʼ����Z�ᴹֱ,��z����
			PixXVec.mult(x, PixXVec.crossProduct(ScreenVec, PixYVec).normalization());	//��ĻX��������Ļ�ᡢ��ĻY��������
			Ray.add(Ray.add(ScreenVec, PixYVec), PixXVec);
			//Paint
			unsigned int color = traceRay(RaySt, Ray, 0);
			g.setPoint(x, y, color);
		}
	}
}
/*--------------------------------[ ׷�ٹ��� ]--------------------------------*/
unsigned int RayTracing::traceRay(Mat<double>& RaySt, Mat<double>& Ray, int level) {
	//Seek Intersecton
	double minDistance = 0, RayFaceTheta, RayFaceThetaTmp;
	Mat<double> intersection;
	Triangle closestTriangle;
	for (int i = 0; i < TriangleSet.size(); i++) {
		double distance = seekIntersection(TriangleSet[i], RaySt, Ray, RayFaceThetaTmp, intersection);
		if (distance > 0 && distance < minDistance) {
			distance < minDistance; closestTriangle = TriangleSet[i]; RayFaceTheta = RayFaceThetaTmp;
		}
	}
	// ���������η��䷽��˥��ϵ�����ٽ��������Ϊ��׼���¼���
	unsigned int color = 0;
	if (level < maxRayLevel) {
		Mat<double> RayTmp;
		color = traceRay(intersection, RayTmp, level + 1);
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
double RayTracing::seekIntersection(Triangle& triangle, Mat<double>& RaySt, Mat<double>& Ray, double& RayFaceTheta, Mat<double>& intersection) {
	//[1]
	Mat<double> edge[2], faceVec, tmp;
	edge[0].add(triangle.p[1], triangle.p[0].negative(edge[0]));
	edge[1].add(triangle.p[2], triangle.p[0].negative(edge[1]));
	faceVec.crossProduct(edge[0], edge[1]);
	//[2]
	double PointFaceDistance = fabs(faceVec.dot(RaySt)) / (faceVec.dot(faceVec));
	RayFaceTheta = faceVec.dot(Ray) / (Ray.norm() * faceVec.norm());
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