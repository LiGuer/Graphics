#ifndef COMPUTATIONAL_GEOMETRY_DELAUNAY_H
#define COMPUTATIONAL_GEOMETRY_DELAUNAY_H

#include <algorithm>
#include <vector>
#include"../Matrix/Matrix.h"

using namespace Matrix;

#define PI 3.141592653589

namespace Geometry {
/*************************************************************************************************
						Delaunay 三角剖分
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
*************************************************************************************************/
Mat<>* Delaunay(Mat<> point[], int n, int& TrianglesNum) {
	std::vector<Mat<>> triAns, triTemp, edgeBuffer;
	std::sort(point, point + n, [](Mat<>& a, Mat<>& b) {				// 将点按坐标x从小到大排序
		return a[0] != b[0] ? a[0] < b[0] : a[1] < b[1];
	});
	//[2]
	Mat<> maxPoint(point[0]), 
				minPoint(point[0]);
	for (int i = 1; i < n; i++) {
		maxPoint = (point[i][0] > maxPoint[0] || (point[i][0] == maxPoint[0] && point[i][1] > maxPoint[1])) ? point[i] : maxPoint;
		minPoint = (point[i][0] < minPoint[0] || (point[i][0] == minPoint[0] && point[i][1] < minPoint[1])) ? point[i] : minPoint;
	}
	Mat<> supertriangle(2, 3), length;
	sub(length, maxPoint, minPoint);
	supertriangle(0, 0) = minPoint[0] - length[0] - 2;    supertriangle(1, 0) = minPoint[1] - 2;
	supertriangle(0, 1) = maxPoint[0] + length[0] + 2;    supertriangle(1, 1) = minPoint[1] - 2;
	supertriangle(0, 2) =(maxPoint[0] + minPoint[0]) / 2; supertriangle(1, 2) = maxPoint[1] + length[1] + 2;
	triTemp.push_back(supertriangle);
	//[3]
	for (int i = 0; i < n; i++) {
		edgeBuffer.clear();
		//[3.2]
		for (int j = 0; j < triTemp.size(); j++) {
			//[3.2.1] 
			Mat<> center, triEdge[3], temp;
			for (int k = 0; k < 3; k++)
				getCol(triTemp[j], k, triEdge[k]);
			double R;
			ThreePoints2Circle(triEdge, center, R);
			double distance = (sub(temp, point[i], center)).norm();
			//[3.2.2]
			if (point[i][0] > center[0] + R) {
				triAns.push_back(triTemp[j]);
				triTemp.erase(triTemp.begin() + j--);
			}
			//[3.2.4]
			else if (distance < R) {
				Mat<> edge(2, 2), p1, p2;
				for (int k = 0; k < 3; k++) {
					getCol(triTemp[j], k, p1);
					getCol(triTemp[j], k + 1) % 3, p2);
					if (p1[0] < p2[0] || (p1[0] == p2[0] && p1[1] < p2[1])) { 
						setCol(edge, 0, p1); 
						setCol(edge, 1, p2);
					}
					else { 
						setCol(edge, 0, p2);
						setCol(edge, 1, p1);
					}
					edgeBuffer.push_back(edge);
				}
				triTemp.erase(triTemp.begin() + j--);
			}
		}
		//[3.3] 
		std::sort(edgeBuffer.begin(), edgeBuffer.end(), [](Mat<> a, Mat<> b) {
			if (a(0, 0) < b(0, 0) || (a(0, 0) == b(0, 0) && a(1, 0) < b(1, 0)))return true;
			if (a(0, 1) < b(0, 1) || (a(0, 1) == b(0, 1) && a(1, 1) < b(1, 1)))return true;
			return false;
		});
		for (int j = 0; j < edgeBuffer.size() - 1; j++) {
			bool flag = 0;
			while (j + 1 < edgeBuffer.size() && edgeBuffer[j] == edgeBuffer[j + 1]) {
				edgeBuffer.erase(edgeBuffer.begin() + j + 1); flag = 1;
			}
			if(flag) { edgeBuffer.erase(edgeBuffer.begin() + j); j--; }
		}
		//[3.4] 
		for (int j = 0; j < edgeBuffer.size(); j++) {
			Mat<> t(2, 3), temp;
			t.setCol(0, edgeBuffer[j].getCol(0, temp)); 
			t.setCol(1, edgeBuffer[j].getCol(1, temp)); 
			t.setCol(2, point[i]);
			triTemp.push_back(t);
		}
	}
	//[4]
	for (int i = 0; i < triTemp.size(); i++) triAns.push_back(triTemp[i]);
	for (int i = 0; i < triAns. size(); i++) {
		Mat<> t;
		for (int j = 0; j < 3; j++) {
			triAns[i].getCol(j, t);
			if (t[0]< minPoint[0] || t[1] < minPoint[1] || t[0] > maxPoint[0] || t[1] > maxPoint[1]) {
				triAns.erase(triAns.begin() + i--); break;
			}
		}
	}
	// [Output]
	TrianglesNum = triAns.size();
	Mat<>* Triangles = (Mat<>*)calloc(TrianglesNum, sizeof(Mat<>));
	for (int i = 0; i < TrianglesNum; i++) Triangles[i] = triAns[i];
	return Triangles;
}


}