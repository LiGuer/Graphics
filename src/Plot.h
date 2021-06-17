/*
Copyright 2020,2021 LiGuer. All Rights Reserved.
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
#ifndef PLOT_H
#define PLOT_H
#include "GraphicsND.h"
//��ѧ��ͼ��
class Plot : public GraphicsND {
public:
	Mat<> pmin, pmax, pdiff, p2v;
	Plot() { 
		clear(0xFFFFFF); g.PaintColor = 0; FACE = 0; LINE = 1;
	};
	void init		(Mat<>& x, Mat<>& y);
	/*---------------- ��ͼ ----------------*/
	void plot		(Mat<>& x, Mat<>& y);
	void plot		(Mat<>& x, Mat<>& y, Mat<>& z);
	void stairs		(Mat<>& y);
	void stairs		(Mat<>& x, Mat<>& y);
	void loglog		(Mat<>& x, Mat<>& y);
	void semilogx	(Mat<>& x, Mat<>& y);
	void semilogy	(Mat<>& x, Mat<>& y);
	void polarplot	(Mat<>& x, Mat<>& y);
	void polarscatter(Mat<>& x, Mat<>& y);
	void contour3	(Mat<>& x);
	void image		(const char* fileName);
	/*---------------- ��ͼ ----------------*/
	void surf		(Mat<>& x, Mat<>& z);
	void mesh		(Mat<>& x, Mat<>& z);
	/*---------------- ��ά���ӻ� ----------------*/
	void coneplot	(Mat<>& x, Mat<>& z);
	void quiver		(Mat<>& x, Mat<>& y);
	void quiver		(Mat<>& x, Mat<>& y, Mat<>& z);
	/*---------------- ���ݷֲ�ͼ ----------------*/
	void histogram	(Mat<>& x);
	void pie		(Mat<>& x, bool* explode = NULL);
	void scatter	(Mat<>& x, Mat<>& y);
	void scatter	(Mat<>& x, Mat<>& y, Mat<>& z);
	void heatmap	(Mat<>& x);
	/*---------------- ��ɢͼ ----------------*/
	void bar		(Mat<>& x);
	void barh		(Mat<>& x);
	void stem		(Mat<>& x);
	/*---------------- ���� ----------------*/
	void show();
	/*---------------- ��� ----------------*/
	void axis		();
	void grid		();
	void title		(const char* words);
};
#endif // !PLOT_H
