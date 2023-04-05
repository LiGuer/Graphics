#ifndef COMPUTATIONAL_GEOMETRY_DELAUNAY_H
#define COMPUTATIONAL_GEOMETRY_DELAUNAY_H

#include <algorithm>
#include <vector>
#include"../Matrix/Matrix.h"

using namespace Matrix;

#define PI 3.141592653589

namespace Geometry {
/*************************************************************************************************
						Delaunay �����ʷ�
*	[����]:
		[1] Delaunay�����ʷ�: ÿ�������ε����Բ�ڲ�����V���κε�
	[����]:
		[1] ���㰴����x��С��������
		[2] ȷ������������
			�����������α�����δȷ���������б� trianglesTemp
		[3] ����ÿһ����
			[3.1] ��ʼ���߻������� edgeBuffer
			[3.2] ���� trianglesTemp �е�ÿһ��������
				[3.2.1] ����������ε�Բ�ĺͰ뾶
				[3.2.2] ����õ������Բ���Ҳ�
					���������ΪDelaunay�����Σ����浽triangles,����temp��ȥ����,����
				[3.2.3] ����õ������Բ�⣨��Ҳ�������Բ�Ҳࣩ
					���������Ϊ��ȷ��,����
				[3.2.4] ����õ������Բ��
					��������β�ΪDelaunay������,�����߱�����edgeBuffer,��temp��ȥ������������
			[3.3] ��edgeBuffer����ȥ��
			[3.4] ��edgeBuffer�еı��뵱ǰ�ĵ������ϳ����������β�������temp triangles��
		[4] ��triangles��trianglesTemp���кϲ�, ����ȥ�볬���������йص�������
*	[Referance]:
		[1] http://paulbourke.net/papers/triangulate/
*************************************************************************************************/
Mat<>* Delaunay(Mat<> point[], int n, int& TrianglesNum) {
	std::vector<Mat<>> triAns, triTemp, edgeBuffer;
	std::sort(point, point + n, [](Mat<>& a, Mat<>& b) {				// ���㰴����x��С��������
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