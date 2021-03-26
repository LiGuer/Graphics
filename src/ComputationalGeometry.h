/*
Copyright 2020 LiGuer. All Rights Reserved.
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
	http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Reference.
[1]Introduction Algorithms.THOMAS H.CORMEN,CHARLES E.LEISERSON,RONALD L.RIVEST,CLIFFORD STEIN
==============================================================================*/
#ifndef COMPUTATIONAL_GEOMETRY_H
#define COMPUTATIONAL_GEOMETRY_H
#include <stdlib.h>
#include <float.h>
#include <algorithm>
#include <vector>
#include <stack>
#include"../LiGu_AlgorithmLib/Mat.h"
namespace ComputationalGeometry{

const double PI = 3.141592653589;
/*----------------[ 2D 二维 ]----------------*/
void ThreePointsToCircle(Mat<double> Points[], Mat<double>& center, double& R);	//平面三点确定圆方程
bool isInCircle(Mat<double> Points[]);											//判断四点共圆
Mat<double>* ConvexHull(Mat<double> point[], int n, int& ansPointNum);			//凸包
Mat<double>* Delaunay(Mat<double> point[], int n, int& TrianglesNum);			//Delaunay三角剖分
/*----------------[ 3D 三维 ]----------------*/
Mat<double>* getSphereFibonacciPoint(int& n);									//球面均匀点分布
/******************************************************************************
*
*                    2D 二维
*
******************************************************************************/
/*--------------------------------[ Segments Intersect 线段相交判断 ]--------------------------------
*	判定条件：
	1.Each segment straddles the line containing the other.
	2.An endpoint of one segment line on the other segment. (the boundary case.)
**---------------------------------------------------------*/
bool isSegmentsIntersect(Mat<double>& a, Mat<double>& b){/*
	double dir_a1 = CrossProduct(a.p[0], a.p[1], b.p[0]),
		dir_a2 = CrossProduct(a.p[0], a.p[1], b.p[1]),
		dir_b1 = CrossProduct(b.p[0], b.p[1], a.p[0]),
		dir_b2 = CrossProduct(b.p[0], b.p[1], a.p[1]);
	if (dir_a1 == 0)
		if (OnSegments_judge(a, b.p[0])) return true; else {}
	else if (dir_a2 == 0)
		if (OnSegments_judge(a, b.p[1])) return true; else {}
	else if (dir_b1 == 0)
		if (OnSegments_judge(b, a.p[0])) return true; else {}
	else if (dir_b2 == 0)
		if (OnSegments_judge(b, a.p[1])) return true; else {}
	else if (dir_a1 != dir_a2 && dir_b1 != dir_b2) return true;
	return false;*/
}
/*--------------------------------[ CircumCircle 三角形外接圆 ]--------------------------------
*	外接圆圆心: 即. 三点确定圆方程问题， 也是任意两边的垂直平分线的交点.直接用 ThreePointsToCircle()方法
**--------------------------------------------------------------------------------------------*/

/*--------------------------------[ ThreePointsToCircle 平面三点确定圆方程 ]--------------------------------
*	[公式]: 圆方程: (x - cx)² + (y - cy)² = R²
*	[算法]: 三点确定圆方程: 即 解行列式:
			| x²+y²    x   y   1 |  =  0
			| x1²+y1²  x1  y1  1 |
			| x2²+y2²  x2  y2  1 |
			| x3²+y3²  x3  y3  1 |
		即.目标三点和圆上(x,y)应该满足方程组:
			(x²+y²)·a + x·b + y·c + 1·d = 0
*	[推导]:
			M11(x²+y²) - M12 x + M13 y - M14 = (x²+y²)·a + x·b + y·c + 1·d = 0
			=> (x² + b/a x) + (y² + c/a y) = - d/a
			=> (x + b/2a)² + (y + c/2a)² = -d/a + b²/4a² + c²/4a²
**------------------------------------------------------------------------------------------------*/
void ThreePointsToCircle(Mat<double> Points[], Mat<double>& center, double& R) {
	center.zero(2, 1);
	Mat<double> A(4, 4);
	for (int i = 0; i < 4; i++)A(i, 3) = 1;
	for (int i = 0; i < 3; i++) {
		double norm = Points[i].norm();
		A(i + 1, 0) = norm * norm;
		A(i + 1, 1) = Points[i][0];
		A(i + 1, 2) = Points[i][1];
	}
	double a = A.comi(0, 0), b = -A.comi(0, 1), c = A.comi(0, 2), d = -A.comi(0, 3);
	center[0] = -b / (2 * a); center[1] = -c / (2 * a);
	R = sqrt(-d / a + b * b / (4 * a * a) + c * c / (4 * a * a));
}
/*--------------------------------[ isInCircle 判断四点共圆 ]--------------------------------
*	[输出]: 圆外-1，圆上0，圆内1
*	三点确定圆方程: 即 解行列式:
		| x1²+y1²  x1  y1  1 | ?= 0
		| x2²+y2²  x2  y2  1 |
		| x3²+y3²  x3  y3  1 |
		| x4²+y4²  x4  y4  1 |
*	[几何解释]: 通过把平面点提升到三维的抛物面中，由于抛物面被平面所截的截面为圆形，四点共面即使共圆，也可以用四面体的体积判断是否共圆。
**------------------------------------------------------------------------------------------------*/
bool isInCircle(Mat<double> Points[]) {
	Mat<double> A(4, 4);
	for (int i = 0; i < 4; i++) {
		A(i, 0) = Points[i].norm();
		A(i, 1) = Points[i][0];
		A(i, 2) = Points[i][1];
		A(i, 4) = 1;
	}
	double t = A.abs();
	return t == 0 ? 1 : 0;
}
/*--------------------------------[ ConvexHull 凸包 ]--------------------------------
*	[算法]: Graham 扫描法
*	[时间复杂度]: O(n logn)
*	[流程]:
		[1] 选择y最小点 p0, 若多个则选其中x最小
		[2] sorted by polar angle in counterclockwise order around p0
			(if more than one point has the same angle, remove all but the one that is farthest from p0)
			* 几何可知，排序后 P1 和最后一点 Pn-1 一定是凸包上的点
		[3] P0,P1 入栈S，P2 为当前点
			if n<2, error "Convex Hull is empty"
		[4] 遍历剩余点 P3 -> Pn-1
			[4.1] while the angle formed by points NEXT-TO-TOP(S),TOP(S),and p makes a nonleft turn
					POP(S)
			[4.2] Pi 入栈
		[5] 最后栈中元素，即结果
*	[Referance]:
		[1] Introduction Algorithms.THOMAS H.CORMEN,CHARLES E.LEISERSON,RONALD L.RIVEST,CLIFFORD STEIN
		[2] Thanks for https://www.cnblogs.com/aiguona/p/7232243.html
**----------------------------------------------------------------------------*/
Mat<double>* ConvexHull(Mat<double> point[], int n, int& ansPointNum){
	// [1] 
	int minCur = 0;
	Mat<double> minPoint(point[0]);
	for (int i = 1; i < n; i++)
		if (point[i][1] < minPoint[1] || (point[i][1] == minPoint[1] && point[i][0] < minPoint[0])) {
			minPoint = point[i]; minCur = i;
		}
	point[0].swap(point[minCur]);
	// [2] 
	std::sort(point + 1, point + n, [&minPoint](Mat<double>& a, Mat<double>& b) {
		if (atan2(a[1] - minPoint[1], a[0] - minPoint[0]) != atan2(b[1] - minPoint[1], b[0] - minPoint[0]))
			return (atan2(a[1] - minPoint[1], a[0] - minPoint[0])) < (atan2(b[1] - minPoint[1], b[0] - minPoint[0]));
		return (a[0] - minPoint[0])* (a[0] - minPoint[0])+ (a[1] - minPoint[1]) * (a[1] - minPoint[1]) 
			< (b[0] - minPoint[0]) * (b[0] - minPoint[0]) + (b[1] - minPoint[1]) * (b[1] - minPoint[1]);
	});
	// [3]
	std::stack<Mat<double>> ConvexHullPoint;
	for (int i = 0; i <= 2; i++)ConvexHullPoint.push(point[i]);
	// [4]
	for (int i = 3; i < n; i++) {
		while (true) {
			Mat<double> prePoint = ConvexHullPoint.top();
			ConvexHullPoint.pop();
			Mat<double> prePointNext = ConvexHullPoint.top(); 
			ConvexHullPoint.push(prePoint);
			// 叉乘判断角度转向
			Mat<double> a, b;
			a.add(prePointNext, prePoint.negative(a));
			b.add(point[i], prePoint.negative(b));
			if (a[0] * b[1] - a[1] * b[0] < 0) break;
			ConvexHullPoint.pop();
		}
		ConvexHullPoint.push(point[i]);
	}
	// [5] Output
	ansPointNum = ConvexHullPoint.size();
	Mat<double>* outputPoint = (Mat<double>*)malloc(ansPointNum * sizeof(Mat<double>));
	memset(outputPoint, 0, ansPointNum * sizeof(Mat<double>));
	for (int i = 0; i < ansPointNum; i++) {
		outputPoint[i] = ConvexHullPoint.top();
		ConvexHullPoint.pop();
	}
	return outputPoint;
}
/*int main() {
	Plot p;
	p.setAxisRange(0, 0, 100, 100);
	Mat<double> point[100],* ConvexHullPoint;
	p.g->PaintColor = 0xFF0000; p.g->PaintSize = 3;
	for (int i = 0; i < 30; i++) {
		point[i].rands(2, 1, 10, 90);
		p.plotPoint(point[i][0], point[i][1]);
	}
	int ConvexHullnum;
	ConvexHullPoint = ConvexHull(point, 30, ConvexHullnum);
	p.g->PaintColor = 0xFFFFFF; p.g->PaintSize = 0;
	for (int i = 0; i < ConvexHullnum; i++) {
		p.plotLine(ConvexHullPoint[i][0], ConvexHullPoint[i][1]
			, ConvexHullPoint[(i + 1) % ConvexHullnum][0], ConvexHullPoint[(i + 1) % ConvexHullnum][1]);
	}p.g->PicWrite("D:/LIGU.ppm");
}*/
/*--------------------------------[ Delaunay 三角剖分 ]--------------------------------
*	[定义]:
		[1] Delaunay三角剖分: 每个三角形的外接圆内不包含V中任何点
	[流程]:
		[1] 将点按坐标x从小到大排序
		[2] 确定超级三角形
			将超级三角形保存至未确定三角形列表 trianglesTemp
		[3] 遍历每一个点
			[3.1] 初始化边缓存数组 edgeBuffer
			[3.2] 遍历 trianglesTemp 中的每一个三角形
				[3.2.1] 计算该三角形的圆心和半径
				[3.2.2] 如果该点在外接圆的右侧
					则该三角形为Delaunay三角形，保存到triangles,并在temp里去除掉,跳过
				[3.2.3] 如果该点在外接圆外（即也不是外接圆右侧）
					则该三角形为不确定,跳过
				[3.2.4] 如果该点在外接圆内
					则该三角形不为Delaunay三角形,将三边保存至edgeBuffer,在temp中去除掉该三角形
			[3.3] 对edgeBuffer进行去重
			[3.4] 将edgeBuffer中的边与当前的点进行组合成若干三角形并保存至temp triangles中
		[4] 将triangles与trianglesTemp进行合并, 并除去与超级三角形有关的三角形
*	[Referance]:
		[1] http://paulbourke.net/papers/triangulate/
**----------------------------------------------------------------------*/
Mat<double>* Delaunay(Mat<double> point[], int n, int& TrianglesNum) {
	std::vector<Mat<double>> triAns, triTemp, edgeBuffer;
	std::sort(point, point + n, [](Mat<double> a, Mat<double> b) {				// 将点按坐标x从小到大排序
		if(a[0] != b[0])return a[0] < b[0]; return a[1] < b[1];
	});
	//[2] 确定超级三角形
	Mat<double> maxPoint(point[0]), minPoint(point[0]);
	for (int i = 1; i < n; i++) {
		maxPoint = (point[i][0] > maxPoint[0] || (point[i][0] == maxPoint[0] && point[i][1] > maxPoint[1])) ? point[i] : maxPoint;
		minPoint = (point[i][0] < minPoint[0] || (point[i][0] == minPoint[0] && point[i][1] < minPoint[1])) ? point[i] : minPoint;
	}
	Mat<double> supertriangle(2, 3), length;
	length.add(maxPoint, minPoint.negative(length));
	supertriangle(0, 0) = minPoint[0] - length[0] - 2; supertriangle(1, 0) = minPoint[1] - 2;
	supertriangle(0, 1) = maxPoint[0] + length[0] + 2; supertriangle(1, 1) = minPoint[1] - 2;
	supertriangle(0, 2) = (maxPoint[0] + minPoint[0]) / 2; supertriangle(1, 2) = maxPoint[1] + length[1] + 2;
	triTemp.push_back(supertriangle);
	//[3] 遍历每一个点
	for (int i = 0; i < n; i++) {
		edgeBuffer.clear();
		//[3.2] 遍历 trianglesTemp 中的每一个三角形
		for (int j = 0; j < triTemp.size(); j++) {
			//[3.2.1] 计算该三角形的圆心和半径
			Mat<double> center, triEdge[3], temp;
			for (int k = 0; k < 3; k++)
				triTemp[j].getCol(k, triEdge[k]);
			double R;
			ThreePointsToCircle(triEdge, center, R);
			double distance = (temp.add(point[i], center.negative(temp))).norm();
			//[3.2.2] 如果该点在外接圆的右侧
			if (point[i][0] > center[0] + R) {
				triAns.push_back(triTemp[j]);
				triTemp.erase(triTemp.begin() + j);
				j--;
			}
			//[3.2.4] 如果该点在外接圆内
			else if (distance < R) {
				Mat<double> edge(2, 2),p1,p2;
				for (int k = 0; k < 3; k++) {
					triTemp[j].getCol(k, p1); triTemp[j].getCol((k + 1) % 3, p2);
					if (p1[0] < p2[0] || (p1[0] == p2[0] && p1[1] < p2[1])) { edge.setCol(0, p1); edge.setCol(1, p2); }
					else { edge.setCol(0, p2); edge.setCol(1, p1); }
					edgeBuffer.push_back(edge);
				}
				triTemp.erase(triTemp.begin() + j);
				j--;
			}
		}
		//[3.3] 对edgeBuffer去重
		std::sort(edgeBuffer.begin(), edgeBuffer.end(),[](Mat<double> a, Mat<double> b) {
			if (a(0, 0) < b(0, 0) || (a(0, 0) == b(0, 0) && a(1, 0) < b(1, 0)))return true;
			if (a(0, 1) < b(0, 1) || (a(0, 1) == b(0, 1) && a(1, 1) < b(1, 1)))return true;
			return false;
		});
		for (int j = 0; j < edgeBuffer.size() - 1; j++) {
			bool flag = 0;
			while (j + 1 < edgeBuffer.size() && edgeBuffer[j] == edgeBuffer[j + 1]) {
				edgeBuffer.erase(edgeBuffer.begin() + j + 1);flag = 1;
			}
			if (flag) { edgeBuffer.erase(edgeBuffer.begin() + j); j--; }
		}
		//[3.4] 将edge buffer中的边与当前的点进行组合成若干三角形并保存至temp triangles中
		for (int j = 0; j < edgeBuffer.size(); j++) {
			Mat<double> t(2, 3), temp;
			t.setCol(0, edgeBuffer[j].getCol(0, temp)); t.setCol(1, edgeBuffer[j].getCol(1, temp)); t.setCol(2, point[i]);
			triTemp.push_back(t);
		}
	}
	//[4] 将triangles与temp triangles进行合并, 并除去与超级三角形有关的三角形
	for (int i = 0; i < triTemp.size(); i++) triAns.push_back(triTemp[i]);
	for (int i = 0; i < triAns.size(); i++) {
		Mat<double> t;
		for (int j = 0; j < 3; j++) {
			triAns[i].getCol(j, t);
			if (t[0]< minPoint[0] || t[1] < minPoint[1] || t[0] > maxPoint[0] || t[1] > maxPoint[1]) {
				triAns.erase(triAns.begin() + i); i--; break;
			}
		}
	}
	// [Output]
	TrianglesNum = triAns.size();
	Mat<double>* Triangles = (Mat<double>*)malloc(TrianglesNum * sizeof(Mat<double>));
	memset(Triangles, 0, TrianglesNum * sizeof(Mat<double>));
	for (int i = 0; i < TrianglesNum; i++) Triangles[i] = triAns[i];
	return Triangles;
}
/*------ example ------
int main() {
	Plot p;
	p.setAxisRange(0, 0, 100, 100);
	Mat<double> point[100],*triangles;
	p.g->PaintColor = 0xFF0000; p.g->PaintSize = 5;
	for (int i = 0; i < 100; i++) {
		point[i].rands(2, 1, 10, 90);
		p.plotPoint(point[i][0], point[i][1]);
	}
	int trianglesnum;
	triangles = Delaunay(point, 100, trianglesnum);
	double x[3], y[3];
	p.g->PaintColor = 0xFFFFFF; p.g->PaintSize = 0;
	for (int i = 0; i < trianglesnum; i++) {
		for (int j = 0; j < 3; j++) {
			x[j] = triangles[i](0, j);
			y[j] = triangles[i](1, j);
		}
		p.ploPolygon(x, y, 3);
	}p.g->PicWrite("D:/LIGU.ppm");
}*/
/******************************************************************************
*
*                    3D 三维
*
******************************************************************************/
/*--------------------------------[ getSphereFibonacciPoint 球面均匀点分布 ]--------------------------------
*	[Referance]:
		[1] Thanks and copyright for https://github.com/SebLague/Boids
**---------------------------------------------------------------------------------------------------------*/
Mat<double>* getSphereFibonacciPoint(int& n) {
	// 球面均匀点分布法
	n = 300;
	Mat<double>* point = (Mat<double>*)malloc(n * sizeof(Mat<double>));
	memset(point, 0, n * sizeof(Mat<double>));
	// 均匀球面点
	double goldenRatio = (1 + sqrt(5)) / 2;				// 黄金分割点
	double angleIncrement = PI * 2 * goldenRatio;
	for (int i = 0; i < n; i++) {
		double t = (double)i / n, inclination = acos(1 - 2 * t), azimuth = angleIncrement * i;
		point[i].zero(3, 1);
		point[i][0] = sin(inclination) * cos(azimuth);
		point[i][1] = sin(inclination) * sin(azimuth);
		point[i][2] = cos(inclination);
	}
	return point;
}
}
#endif