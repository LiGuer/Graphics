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
#ifndef GRAPHICS_FILECODE_H
#define GRAPHICS_FILECODE_H
#include "../../LiGu_AlgorithmLib/Mat.h"
#include "RGB.h"
namespace GraphicsFileCode {
/******************************************************************************
*					.PPM 文件编码/解码
*	[格式]:
		[1] P6 + 图片格式 + 宽高 + 最大像素值
		[2] RGB像素数据
******************************************************************************/
void ppmWrite(const char* fileName, Mat<RGB>& image) {
	FILE* fp = fopen(fileName, "wb");
	fprintf(fp, "P6\n%d %d\n255\n", image.cols, image.rows);			// 写图片格式、宽高、最大像素值
	fwrite(image.data, 1, image.size() * 3, fp);						// 写RGB数据
	fclose(fp);
}
/******************************************************************************
*					.STL 文件编码/解码
******************************************************************************/
void stlFileRead(const char* fileName, Mat<>& faceVec, Mat<>& p1, Mat<>& p2, Mat<>& p3) {
	FILE* fin = fopen(fileName, "rb");
	unsigned char head[80];
	unsigned int  N; 
	short t;
	float p[12];
	fread(head, 80, 1, fin);
	fread(&N,    4, 1, fin);
	faceVec.zero(3, N);
	p1.		zero(3, N);
	p2.		zero(3, N);
	p3.		zero(3, N);
	for (int i = 0; i < N; i++) {
		fread(p, 12 * 4, 1, fin);
		fread(&t, 2    , 1, fin);
		faceVec	(0, i) = p[0];	faceVec	(1, i) = p[1];	faceVec	(2, i) = p[2];
		p1		(0, i) = p[3];	p1		(1, i) = p[4];	p1		(2, i) = p[5];
		p2		(0, i) = p[6];	p2		(1, i) = p[7];	p2		(2, i) = p[8];
		p3		(0, i) = p[9];	p3		(1, i) = p[10];	p3		(2, i) = p[11];
	}
}
}
#endif

