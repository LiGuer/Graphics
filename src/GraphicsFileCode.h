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
static void ppmRead(const char* fileName, Mat<RGB>& image) {
	FILE* fi = fopen(fileName, "rb");
	int rows, cols;
	fscanf(fi, "P6\n%d %d\n255\n", &cols, &rows);						// 读图片格式、宽高、最大像素值
	image.alloc(rows, cols);
	fread(image.data, 1, image.size() * 3, fi);							// 读RGB数据
	fclose(fi);
}
static void ppmWrite(const char* fileName, Mat<RGB>& image) {
	FILE* fo = fopen(fileName, "wb");
	fprintf(fo, "P6\n%d %d\n255\n", image.cols, image.rows);			// 写图片格式、宽高、最大像素值
	fwrite(image.data, 1, image.size() * 3, fo);						// 写RGB数据
	fclose(fo);
}
/******************************************************************************
*					.STL 文件编码/解码
*	[格式]:
		[1] 头			(80B)
		[2] 三角面个数	(4B)
		[3] 对每个三角面(50B)(循环)
			[3.1] 法向量(3x4B)
			[3.2] 顶点1 (3x4B)
			[3.3] 顶点2 (3x4B)
			[3.4] 顶点3 (3x4B)
			[3.5] 属性	(2B)
******************************************************************************/
static void stlRead(const char* fileName, Mat<>& faceVec, Mat<>& p1, Mat<>& p2, Mat<>& p3, Mat<short>& attribute) {
	FILE* fi = fopen(fileName, "rb");
	unsigned char head[80];
	unsigned int  N; 
	float p[12];
	fread(head, 80, 1, fi);
	fread(&N,    4, 1, fi);
	faceVec.alloc(3, N);
	p1.		alloc(3, N);
	p2.		alloc(3, N);
	p3.		alloc(3, N);
	attribute.alloc(N);
	for (int i = 0; i < N; i++) {
		fread(p, 12 * 4, 1, fi);
		fread(attribute.data + i, 2, 1, fi);
		faceVec	(0, i) = p[0];	faceVec	(1, i) = p[1];	faceVec	(2, i) = p[2];
		p1		(0, i) = p[3];	p1		(1, i) = p[4];	p1		(2, i) = p[5];
		p2		(0, i) = p[6];	p2		(1, i) = p[7];	p2		(2, i) = p[8];
		p3		(0, i) = p[9];	p3		(1, i) = p[10];	p3		(2, i) = p[11];
	}
	fclose(fi);
}
static void stlWrite(const char* fileName, const char* head, Mat<float>& faceVec, Mat<float>& p1, Mat<float>& p2, Mat<float>& p3, Mat<short>& attribute) {
	FILE* fo = fopen(fileName, "wb");
	unsigned int  N = p1.cols;
	fwrite(head, 80, 1, fo);
	fwrite(&N,    4, 1, fo);
	for (int i = 0; i < N; i++) {
		fwrite(&faceVec(0, i), 4, 1, fo); fwrite(&faceVec(1, i), 4, 1, fo); fwrite(&faceVec(2, i), 4, 1, fo);
		fwrite(&p1	   (0, i), 4, 1, fo); fwrite(&p1	 (1, i), 4, 1, fo); fwrite(&p1	   (2, i), 4, 1, fo);
		fwrite(&p2	   (0, i), 4, 1, fo); fwrite(&p2	 (1, i), 4, 1, fo); fwrite(&p2	   (2, i), 4, 1, fo);
		fwrite(&p3	   (0, i), 4, 1, fo); fwrite(&p3	 (1, i), 4, 1, fo); fwrite(&p3	   (2, i), 4, 1, fo);
		fwrite(attribute.data + i, 2, 1, fo);
	}
	fclose(fo);
}

}
#endif

