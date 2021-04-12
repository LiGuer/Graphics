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
==============================================================================*/
#ifndef RAY_TRACING_H
#define RAY_TRACING_H
#include "Graphics.h"
#define PI 3.141592653589

class RayTracing {
public:
	struct RGB { 
		unsigned char R, G, B; 
		RGB(unsigned int a) { *this = a; }
		RGB& operator=(const RGB& a) { R = a.R; G = a.G; B = a.B; return *this; }
		RGB& operator=(unsigned int& a) { 
			R = a >> 16; 
			G = a >> 8; 
			B = a;
			return *this;
		}
		RGB& operator*=(double a) {
			R *= a; G *= a; B *= a;
			return *this;
		}
	};
	struct Material {															//����
		RGB color = 0;
		double reflectance = 1, refractiveIndex = 0;
	};
	struct Triangle {															//������
		Mat<double> p[3];
		Material* material = NULL;
	};
	/*---------------- �������� ----------------*/
	Graphics g;																	//����ͼ��ѧ��
	Mat<double> Eye{ 3,1 }, gCenter{ 3,1 };
	int maxRayLevel = 20;
	double refractiveIndexBuffer = 1;
	std::vector<Triangle> TriangleSet;											//�����μ�
	std::vector<Material> MaterialSet;											//���ʼ�
	std::vector<Mat<double>> LightSource;										//���ʼ�
	/*---------------- �ײ� ----------------*/
	RayTracing() { ; }
	RayTracing(int width, int height) { init(width, height); }					//���캯��
	~RayTracing() { ; }															//��������
	void init(int width, int height);											//��ʼ��
	/*---------------- DRAW ----------------*/
	void paint();																		//��Ⱦ
	RGB traceRay(Mat<double>& RaySt, Mat<double>& Ray, RGB& color, int level);
	double seekIntersection(Triangle& triangle, Mat<double>& RaySt, Mat<double>& Ray, Mat<double>& FaceVec, Mat<double>& intersection);	//�󽻵�
	//2-D
	void drawTriangle(Mat<double>& p1, Mat<double>& p2, Mat<double>& p3, Material* material = NULL);							//��������
	void drawRectangle(Mat<double>& sp, Mat<double>& ep, Mat<double>* direct = NULL);											//������
	void drawQuadrilateral(Mat<double>& p1, Mat<double>& p2, Mat<double>& p3, Mat<double>& p4, Material* material = NULL);		//���ı���
	void drawPolygon(Mat<double> p[], int n, Material* material = NULL);														//�������
	void drawCircle(Mat<double>& center, double r, Mat<double>* direct = NULL);						//��Բ
	void drawEllipse(Mat<double>& center, double rx, double ry, Mat<double>* direct = NULL);		//����Բ
	void drawSurface(Mat<double> z, double xs, double xe, double ys, double ye);					//������	
	// 3-D
	void drawTetrahedron(Mat<double>& p1, Mat<double>& p2, Mat<double>& p3, Mat<double>& p4, Material* material = NULL);		//��������
	void drawCuboid(Mat<double>& pMin, Mat<double>& pMax, Material* material = NULL);											//������
	void drawPolyhedron(Mat<double>* p, int n, Material* material = NULL);														//��������
	void drawFrustum(Mat<double>& st, Mat<double>& ed, double Rst, double Red, double delta = 5, Material* material = NULL);	//��Բ̨
	void drawCylinder(Mat<double>& st, Mat<double>& ed, double r, double delta = 5, Material* material = NULL);					//��Բ��
	void drawSphere(Mat<double>& center, double r, Material* material = NULL);													//����
	void drawSphere2(Mat<double>& center, double r, int n = 300, Material* material = NULL);									//����
	void drawEllipsoid(Mat<double>& center, Mat<double>& r, Material* material = NULL);											//������
	void drawBody(Mat<double>& center, Mat<double>& r, Material* material = NULL);												//������
	void drawBezierBody(Mat<double> p[], int n, Material* material = NULL);
};

#endif