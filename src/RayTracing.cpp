#include "RayTracing.h"
/*--------------------------------[ 初始化 ]--------------------------------*/
void RayTracing::init(int width, int height) {
	g.init(width, height);
}
/*--------------------------------[ 渲染 ]--------------------------------*/
void RayTracing::paint() {
	double minDistance = 0, RayFaceDistance;
	
	Mat<double> ScreenVec, PixYVec(3, 1), PixXVec, Ray, RaySt;
	ScreenVec.add(gCenter, ScreenVec.negative(Eye));
	for (int x = 0; x < g.Canvas.rows; x++) {
		for (int y = 0; y < g.Canvas.cols; y++) {
			//Compute Ray
			double t[] = { 1,-ScreenVec[0] / ScreenVec[1],0 };
			PixYVec.mult(y, PixYVec.getData(t));								//屏幕Y向轴始终与Z轴垂直,无z分量
			PixXVec.mult(x, PixXVec.crossProduct(ScreenVec, PixYVec).normalization());	//屏幕X向轴与屏幕轴、屏幕Y向轴正交
			Ray.add(Ray.add(ScreenVec, PixYVec), PixXVec);
			//Paint
			unsigned int color = traceRay(RaySt, Ray, 0);
			g.setPoint(x, y, color);
		}
	}
}
/*--------------------------------[ 追踪光线 ]--------------------------------*/
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
	// 计算三角形反射方向及衰减系数，再将反射光线为基准从新计算
	unsigned int color = 0;
	if (level < maxRayLevel) {
		Mat<double> RayTmp;
		color = traceRay(intersection, RayTmp, level + 1);
	}
	return color;
}
/*--------------------------------[ 求交点 ]--------------------------------
*	[流程]:
		[1] 计算三角形所在面矢量
		[2] 计算点线距离、光线与面矢量夹角、光线与面相交所走过的距离
		[3] 计算光线与面的交点
		[4] 判断交点是否在三角形内部, 若否返回-1
*	[算法]:
		
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