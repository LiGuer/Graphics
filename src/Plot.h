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
//科学制图类
class Plot {
public:
	GraphicsND G;
	Mat<> pmin, pmax, pdiff, p2v;
	Plot() { 
		G.init(700, 500); G.clear(0xFFFFFF); G.g.PaintColor = 0; G.FACE = 0; G.LINE = 1;
	};
	void init		(Mat<>& x, Mat<>& y);
	/*---------------- 线图 ----------------*/
	void plot		(Mat<>& x, Mat<>& y);
	void plot		(Mat<>& x, Mat<>& y, Mat<>& z);
	void stairs		(Mat<>& y);
	void stairs		(Mat<>& x, Mat<>& y);
	void loglog		(Mat<>& x, Mat<>& y);
	void semilogx	(Mat<>& x, Mat<>& y);
	void semilogy	(Mat<>& x, Mat<>& y);
	void polarplot	(Mat<>& x, Mat<>& y);
	void contour	(Mat<>& x, Mat<>& y);
	/*---------------- 面图 ----------------*/
	void surf		(Mat<>& x, Mat<>& z);
	void mesh		(Mat<>& x, Mat<>& z);
	/*---------------- 三维可视化 ----------------*/
	void coneplot		(Mat<>& x, Mat<>& z);
	/*---------------- 数据分布图 ----------------*/
	void histogram	(Mat<>& x);
	void pie		(Mat<>& x);
	void scatter	(Mat<>& x, Mat<>& y);
	void scatter	(Mat<>& x, Mat<>& y, Mat<>& z);
	void hetMap		(Mat<>& x);
	/*---------------- 离散图 ----------------*/
	void bar		(Mat<>& x);
	void barh		(Mat<>& x);
	void stem		(Mat<>& x);
	/*---------------- 功能 ----------------*/
	void show();
	/*---------------- 标记 ----------------*/
	void axis		();
	void grid		();
	void title		(const char* words);
};
#endif // !PLOT_H
