#ifndef GRAPHICS_IO_H
#define GRAPHICS_IO_H

#include <vector>
#include "../../../Math/src/Matrix/Mat.h"
#include "RGB.h"

namespace GraphicsIO {
/******************************************************************************
*					.PPM 文件编码/解码
*	[格式]:
		[1] P() + 图片格式 + 宽高 + 最大像素值
		[2] RGB像素数据
*	[头格式]:
		P1	Bitmap	ASCII			P4	Bitmap	Binary
		P2	Graymap	ASCII			P5	Graymap	Binary
		P3	Pixmap	ASCII			P6	Pixmap	Binary
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

static void ppmWrite(const char* fileName, Mat<ARGB>& image) {
	Mat<RGB> imgT(image.rows, image.cols);
	for (int i = 0; i < image.size(); i++) {
		imgT(i).B = (unsigned char) image(i);
		imgT(i).G = (unsigned char)(image(i) >> 8);
		imgT(i).R = (unsigned char)(image(i) >>16);
	}

	ppmWrite(fileName, imgT);
}

static void ppmWrite(const char* fileName, Mat<unsigned char>& image) {
	FILE* fo = fopen(fileName, "wb");
	fprintf(fo, "P5\n%d %d\n255\n", image.cols, image.rows);			// 写图片格式、宽高、最大像素值
	fwrite(image.data, 1, image.size(), fo);							// 写RGB数据
	fclose(fo);
}

static void ppmWrite(const char* fileName, Mat<double>& image) {
	Mat<unsigned char> t(image.rows, image.cols);
	for (int i = 0; i < image.size(); i++) 
		t[i] = image[i] * 0xFF;
	ppmWrite(fileName, t);
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
	p1.alloc(3, N);
	p2.alloc(3, N);
	p3.alloc(3, N);
	attribute.alloc(N);

	for (int i = 0; i < N; i++) {
		fread(p, 12 * 4, 1, fi);
		fread(attribute.data + i, 2, 1, fi);
		faceVec	(0, i) = p[0];
		faceVec	(1, i) = p[1];
		faceVec	(2, i) = p[2];

		p1(0, i) = p[3];
		p1(1, i) = p[4];
		p1(2, i) = p[5];
		p2(0, i) = p[6];
		p2(1, i) = p[7];
		p2(2, i) = p[8];
		p3(0, i) = p[9];
		p3(1, i) = p[10];
		p3(2, i) = p[11];
	}
	fclose(fi);
}

static void stlWrite(
	const char* fileName, 
	const char* head,
	Mat<float>& faceVec, 
	Mat<float>& p1, 
	Mat<float>& p2, 
	Mat<float>& p3, 
	Mat<short>& attribute
) {
	FILE* fo = fopen(fileName, "wb");
	unsigned int  N = p1.cols;
	fwrite(head, 80, 1, fo);
	fwrite(&N,    4, 1, fo);

	for (int i = 0; i < N; i++) {
		fwrite(&faceVec(0, i), 4, 1, fo); 
		fwrite(&faceVec(1, i), 4, 1, fo); 
		fwrite(&faceVec(2, i), 4, 1, fo);

		fwrite(&p1(0, i), 4, 1, fo);
		fwrite(&p1(1, i), 4, 1, fo);
		fwrite(&p1(2, i), 4, 1, fo);
		fwrite(&p2(0, i), 4, 1, fo);
		fwrite(&p2(1, i), 4, 1, fo);
		fwrite(&p2(2, i), 4, 1, fo);
		fwrite(&p3(0, i), 4, 1, fo);
		fwrite(&p3(1, i), 4, 1, fo);
		fwrite(&p3(2, i), 4, 1, fo);

		fwrite(attribute.data + i, 2, 1, fo);
	}
	fclose(fo);
}

static void stlWrite(
	const char* fileName,
	vector<vector<double>>& triangleSet
) {
	int n = triangleSet.size();
	char head[80];
	Mat<float>
		faceVec(3, n),
		p1(3, n),
		p2(3, n),
		p3(3, n);
	Mat<short>
		attribute(n);

	for (int i = 0; i < n; i++) {
		p1(0, i) = triangleSet[i][0];
		p1(1, i) = triangleSet[i][1];
		p1(2, i) = triangleSet[i][2];
		p2(0, i) = triangleSet[i][3];
		p2(1, i) = triangleSet[i][4];
		p2(2, i) = triangleSet[i][5];
		p3(0, i) = triangleSet[i][6];
		p3(1, i) = triangleSet[i][7];
		p3(2, i) = triangleSet[i][8];
	}

	stlWrite(fileName, head, faceVec, p1, p2, p3, attribute);
}

}
#endif