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
#ifndef DIGITAL_IMAGE_PROCESSING_H
#define DIGITAL_IMAGE_PROCESSING_H
#include <stdlib.h>
#include <time.h>
#define PI 3.141592653589

namespace DigitalImageProcessing {

struct RGB { double R, G, B; };
/*************************************************************************************************
*								图像输入/输出
*	[定义]: Zn+1 = Zn² + C
			所有能使Zn+1不发散的复数点C, 所构成的集合,即 Mandelbrot Set
			(不发散,不一定收敛,有可能在几个不同点来回跳)
*	[性质]: |Zn|>2不可能收敛, 即Mandelbrot Set在半径为2的圆内.
*************************************************************************************************/
/*
#include "../../LiGu_AlgorithmLib/Mat.h"
#include "../../LiGu_AlgorithmLib/BasicMachineLearning.h"
#include "../opencv2-include/opencv2/opencv.hpp"
#pragma comment(lib,"opencv2-include/opencv_world430.lib")
Mat<>* Input(const char* inputImgUrl, Mat<>* data) {
	cv::Mat input = cv::imread(inputImgUrl, cv::IMREAD_COLOR);
	for (int k = 0; k < 3; k++) data[k].zero(input.rows, input.cols);

	for (int i = 0; i < input.rows * input.cols; i++) {
		cv::MatIterator_<cv::Vec3b> it = input.begin<cv::Vec3b>() + i;
		for (int k = 0; k < 3; k++)  data[2 - k][i] = (double)(*it)[k] / 255;
	}
	return data;
}
void Output(const char* outputImgUrl, Mat<>* data) {
	unsigned char* output = (unsigned char*)calloc(data[0].size() * 3, sizeof(unsigned char));
	for (int i = 0; i < data[0].size(); i++)
		for (int k = 0; k < 3; k++)
			output[i * 3 + k] = (data[k])[i] * 255;
	FILE* fp = fopen(outputImgUrl, "wb");
	fprintf(fp, "P6\n%d %d\n255\n", data[0].cols, data[0].rows);	// 写图片格式、宽高、最大像素值
	fwrite(output, 1, data[0].size() * 3, fp);	// 写RGB数据
	fclose(fp);
	free(output);
}
void Output(const char* outputImgUrl, Mat<>& data) {
	unsigned char* output = (unsigned char*)calloc(data.size(), sizeof(unsigned char));
	for (int i = 0; i < data.size(); i++)
			for (int k = 0; k < 3; k++)
				output[i] = data[i] * 255;
	FILE* fp = fopen(outputImgUrl, "wb");
	fprintf(fp, "P5\n%d %d\n255\n", data.cols, data.rows);	// 写图片格式、宽高、最大像素值
	fwrite(output, 1, data.size(), fp);	// 写RGB数据
	fclose(fp);
	free(output);
}*/
/*************************************************************************************************
*								二值化
*	[目的]: 基于阈值，将图像简化为纯黑白图.
*	[公式]: BinarizationImage = Image > threshold ? 1 : 0;
*************************************************************************************************/
Mat<>& Binarization(Mat<>& input, Mat<>& output, double threshold = 0.5) {
	output.zero(input.rows, input.cols);
	for (int i = 0; i < input.size(); i++)
		output[i] = input[i] > threshold ? 1 : 0;
	return output;
}
/*************************************************************************************************
*								颜色聚类
*	[目的]: 简化聚类图像中的色彩.
*	[算法]: K-Mean均值聚类
*************************************************************************************************/
Mat<>* ColorCluster(Mat<>* input, Mat<>* output, int K = 3, int TimesMax = 0x7FFFFFFF) {
	// Process input & output
	Mat<> data(3, input[0].size());
	for (int k = 0; k < 3; k++)
		for (int i = 0; i < input[0].rows; i++)
			for (int j = 0; j < input[0].cols; j++)
				data(k, i * input[0].cols + j) = (input[k])(i, j);
	for (int k = 0; k < 3; k++) output[k].zero(input[0].rows, input[0].cols);
	// Color Cluster
	time_t now; srand((unsigned)time(&now));
	Mat<> Center;
	Mat<int> Cluster, Cluster_Cur;
	BasicMachineLearning::K_Mean(data, K, Center, Cluster, Cluster_Cur, TimesMax);
	for (int i = 0; i < K; i++)
		for (int j = 0; j < Cluster_Cur[i]; j++)
			for (int dim = 0; dim < 3; dim++)
				(output[dim])(Cluster(i, j) / output[0].cols, Cluster(i, j) % output[0].cols) = Center(dim, i);
	return output;
}
/*************************************************************************************************
*								边缘检测
*	[目的]: 标识数字图像中亮度变化明显的点.
*	[公式]: EdgeImage = Conv(Image , SobelKernel)
*************************************************************************************************/
Mat<>& EdgeDetection(Mat<>& input, Mat<>& output) {
	Mat<> SobelKernel(3, 3);
	{
		double t[] = {
			-1,0,1,
			-2,0,2,
			-1,0,1
		}; SobelKernel.getData(t);
	}
	Mat<> output_x, output_y;
	output_x.conv(input, SobelKernel, 1);
	output_y.conv(input, SobelKernel.transpose(output_y), 1);
	output.zero(input.rows, input.cols);
	for (int i = 0; i < input.size(); i++)
		output[i] = sqrt(output_x[i] * output_x[i] + output_y[i] * output_y[i]);
	return output;
}
/*************************************************************************************************
*								傅里叶变换
*	[目的]: 转频域图像.
*************************************************************************************************/
Mat<>& FourierTransform(Mat<>& input, Mat<>& output) {
	return output;
}
Mat<>& InvFourierTransform(Mat<>& input, Mat<>& output) {
	return output;
}
/*************************************************************************************************
*								Gauss 滤波
* [输入]: input: 输入原图 dst: 模糊图像  size: 核的大小  sigma: 正态分布标准差
*************************************************************************************************/
Mat<>& GaussFilter(Mat<>& input, int size, float sigma, Mat<>& output) {
	if (size <= 0 || sigma == 0)return output;
	//二维Gauss核生成
	Mat<> kernel(size, size);
	double sum = 0;
	for (int y = 0; y < size; y++) {
		for (int x = 0; x < size; x++) {
			kernel(x, y) = (1 / (2 * PI * sigma * sigma))
				* exp(-((x - size / 2) * (x - size / 2) + (y - size / 2) * (y - size / 2)) / (2 * sigma * sigma));
			sum += kernel(x, y);
		}
	}
	kernel.mult(1 / sum, kernel);
	//Gauss卷积
	output.conv(input, kernel, 1);
	return output;
}
/*************************************************************************************************
*								转灰度图
*	[目的]: RGB三通道合并为灰度一通道
*	[公式]: Gray = 0.3 R + 0.59 G + 0.11 B
*************************************************************************************************/
Mat<>& Gray(Mat<>* input, Mat<>& output, double Rk = 0.3, double Gk = 0.59, double Bk = 0.11) {
	output.zero(input[0].rows, input[0].cols);
	Mat<> tmp;
	output.add(output, tmp.mult(Rk / (Rk + Gk + Bk), input[0]));
	output.add(output, tmp.mult(Gk / (Rk + Gk + Bk), input[1]));
	output.add(output, tmp.mult(Bk / (Rk + Gk + Bk), input[2]));
	return output;
}
/*************************************************************************************************
*								直方图
* [目的]: 统计[0,255]亮度的像素个数分布.
*************************************************************************************************/
Mat<int>& Histograms(Mat<>& input, Mat<>& output) {
	output.zero(1, 255);
	for (int i = 0; i < input.rows; i++)
		for (int j = 0; j < input.cols; j++)
			output[(unsigned char)(input(i, j) * 255)]++;
}
Mat<int>& Histograms(Mat<>* input, Mat<>& output) {
	output.zero(3, 255);
	for (int k = 0; k < 3; k++)
		for (int i = 0; i < input->rows; i++)
			for (int j = 0; j < input->cols; j++)
				output(k, (unsigned char)(input[k](i, j) * 255))++;
}
/*--------------------------------[ 反相 ]--------------------------------
[目的]: 所有颜色换成其补色
[公式]: InvImage = 1 - Image
-------------------------------------------------------------------------*/
Mat<>& Invert(Mat<>& input, Mat<>& output) {
	output.mult(-1, input); return output;
}
Mat<>* Invert(Mat<>* input, Mat<>* output) {
	for (int k = 0; k < 3; k++)  Invert(input[k], output[k]); return output;
}
}
/*//Example
int main() {
	Mat<> img[3];
	DigitalImageProcessing::Input("IMG01.jpg", img);
	Mat<> colorclusterImg[3];
	DigitalImageProcessing::ColorCluster(img, colorclusterImg, 3, 20);
	Mat<> grayImg;
	DigitalImageProcessing::Gray(colorclusterImg, grayImg);
	Mat<> edgedetectionImg;
	DigitalImageProcessing::EdgeDetection(grayImg, edgedetectionImg);
	DigitalImageProcessing::Output("img.ppm", edgedetectionImg);
}
*/
#endif