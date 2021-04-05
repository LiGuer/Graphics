#include "RayTracing.h"
/*--------------------------------[ 初始化 ]--------------------------------*/
void RayTracing::init(int width, int height) {
	g.init(width, height);
}
/*--------------------------------[ 渲染 ]--------------------------------
*	[过程]:
		[1] 计算屏幕矢量、屏幕X,Y向轴
		[2] 对屏幕每个像素遍历
			[3] 计算像素矢量、光线矢量、光线追踪起点
			[4] 光线追踪算法
			[5] 基于结果绘制该像素色彩
-------------------------------------------------------------------------*/
void RayTracing::paint() {
	//[1]
	Mat<double> ScreenVec, ScreenXVec, ScreenYVec(3, 1);
	ScreenVec.add(gCenter, Eye.negative(ScreenVec));												//屏幕轴由眼指向屏幕中心
	{ double t[] = { 1,-ScreenVec[0] / ScreenVec[1],0 }; ScreenYVec.getData(t).normalization(); }	//屏幕Y向轴始终与Z轴垂直,无z分量
	ScreenXVec.crossProduct(ScreenVec, ScreenYVec).normalization();									//屏幕X向轴与屏幕轴、屏幕Y向轴正交
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
/*--------------------------------[ 追踪光线 ]--------------------------------
*	[过程]:
		[1] 遍历三角形集合中的每一个三角形
			[2]	判断光线和该三角形是否相交、光线走过距离、交点坐标、光线夹角
			[3] 保留光线走过距离最近的三角形的相关数据
		[4] 如果该光线等级小于设定的阈值等级
			计算三角形反射方向，将反射光线为基准重新计算
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
/*--------------------------------[ 求交点 ]--------------------------------
*	[流程]:
		[1] 计算三角形所在面矢量
		[2] 计算点线距离、光线与面矢量夹角、光线与面相交所走过的距离
		[3] 计算光线与面的交点
		[4] 判断交点是否在三角形内部, 若否返回-1
*	[算法]:
		
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