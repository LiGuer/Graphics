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
#ifndef ARCHITECTURE_H
#define ARCHITECTURE_H
#include "../src/GraphicsND.h"
namespace Architecture {
/******************************************************************************
*                    楼梯
******************************************************************************/
void Stairs_1(Mat<>& zero, double Len, double Height,int num, GraphicsND* G) {
	Mat<> p1(3), p2(3);
	G->drawCuboid(
		p1.getData(-Len / 2,-Len / 2, -0.1) += zero,
		p2.getData( Len / 2,-Len / 4,    0) += zero
	);
	G->drawCuboid(
		p1.getData(-Len / 2, Len / 4, Height / 2 - 0.1) += zero,
		p2.getData( Len / 2, Len / 2, Height / 2      ) += zero
	);
	p1 = zero; p1[0] += Len / 4;
	G->drawStairs(p1, Len / 2, Len / 2, Height / 2, num);
	p1 = zero; p1[0] -= Len/4; p1[2] += Height / 2;
	G->drawStairs(p1,   Len / 2, Len / 2, Height / 2, num);
}
void Stairs_2(Mat<>& zero, double Len, double Height,int num, GraphicsND* G) {
	Mat<> p1(3), p2(3);
	G->drawCuboid(
		p1.getData(-Len / 2,-Len / 2, -0.1) += zero,
		p2.getData( Len / 2,-Len / 4,    0) += zero
	);
	G->drawCuboid(
		p1.getData(-Len / 2, Len / 4, Height / 2 - 0.1) += zero,
		p2.getData( Len / 2, Len / 2, Height / 2      ) += zero
	);
	p1 = zero; 
	G->drawStairs(p1, Len / 2, Len / 2, Height / 2, num);
	p1 = zero; p1[0] -= Len / 8; p1[2] += Height / 2;
	G->drawStairs(p1, Len / 4, Len / 2, Height / 2, num);
	p1 = zero; p1[0] += Len / 8; p1[2] += Height / 2;
	G->drawStairs(p1, Len / 4, Len / 2, Height / 2, num);
}
/******************************************************************************
*                    镂空平面
******************************************************************************/
void Wall(Mat<>& st, Mat<>& ed, Mat<>& holeSt, Mat<>& holeEd, GraphicsND* G) {
	Mat<> p1(3);
	(p1 = ed)[0] = holeSt[0];  G->drawCuboid(st, p1);
	(p1 = ed)[1] = holeSt[1];  G->drawCuboid(st, p1);
	(p1 = st)[0] = holeEd[0];  G->drawCuboid(p1, ed);
	(p1 = st)[1] = holeEd[1];  G->drawCuboid(p1, ed);
}
/******************************************************************************
*                    栏杆
******************************************************************************/
void Handrail(Mat<>& st, Mat<>& ed, double Height, double Delta, GraphicsND* G) {
	Mat<> p1 = st, p2 = ed, deltaP;
	p1[2] += Height;
	p2[2] += Height;
	deltaP.sub(ed, st).normalized() *= Delta;
	G->drawPipe(p1, p2, 0.05);
	p1 = st;
	while (p1[0] <= ed[0] && p1[1] <= ed[1] && p1[2] <= ed[2]) {
		p2 = p1; p2[2] += Height;
		G->drawPipe(p1, p2, 0.05);
		p1 += deltaP;
	}
}
}
#endif