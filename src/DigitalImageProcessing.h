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
#include "../../LiGu_AlgorithmLib/Mat.h"
#include "../../LiGu_AlgorithmLib/BasicMachineLearning.h"
#define PI 3.141592653589
namespace DigitalImageProcessing {
/*************************************************************************************************
*								基础操作
*	二值化	:	基于阈值，将图像简化为纯黑白图	:	BinarizationImage = Image > threshold ? 1 : 0;
*	反相	:	所有颜色换成其补色				:	InvImage = 1 - Image
*	转灰度图:	多通道(RGB)加权合并为灰度一通道	:	Gray = 0.3 R + 0.59 G + 0.11 B
*************************************************************************************************/
Mat<>& Binarization	(Mat<>& in, Mat<>& out, double threshold = 0.5) {
	return out.zero(in).function(in, [&threshold](double x) { x > threshold ? 1.0 : 0.0; });
}
Mat<>& Invert		(Mat<>& in, Mat<>& out) { return out.mul(-1, in); }
Mat<>& Gray			(Mat<>* in, Mat<>& out, double Rk = 0.3, double Gk = 0.59, double Bk = 0.11) {
	out.zero(in[0]);
	Mat<> t;
	out += t.mul(Rk / (Rk + Gk + Bk), in[0]);
	out += t.mul(Gk / (Rk + Gk + Bk), in[1]);
	out += t.mul(Bk / (Rk + Gk + Bk), in[2]);
	return out;
}
Mat<>& Gray(Mat<>* in, Mat<>& out, double* rate, int N) {
	out.zero(in[0]);
	Mat<> t;
	for (int i = 0; i < N; i++) out += t.mul(rate[i], in[i]);
	return out;
}
/*************************************************************************************************
*								颜色聚类
*	[目的]: 简化聚类图像中的色彩.
*	[算法]: K-Mean均值聚类
*************************************************************************************************/
Mat<>* ColorCluster(Mat<>* in, Mat<>* out, int K = 3, int TimesMax = 0x7FFFFFFF) {
	// Process in & out
	Mat<> data(3, in[0].size());
	for (int k = 0; k < 3; k++)
		for (int i = 0; i < in[0].rows; i++)
			for (int j = 0; j < in[0].cols; j++)
				data(k, i * in[0].cols + j) = (in[k])(i, j);
	for (int k = 0; k < 3; k++) out[k].zero(in[0].rows, in[0].cols);
	// Color Cluster
	Mat<> Center;
	Mat<int> Cluster, Cluster_Cur;
	//BasicMachineLearning::K_Mean(data, K, Center, Cluster, Cluster_Cur, TimesMax);
	for (int i = 0; i < K; i++)
		for (int j = 0; j < Cluster_Cur[i]; j++)
			for (int dim = 0; dim < 3; dim++)
				(out[dim])(Cluster(i, j) / out[0].cols, Cluster(i, j) % out[0].cols) = Center(dim, i);
	return out;
}
/*************************************************************************************************
*								特殊算子
*************************************************************************************************/
double SobelKernelTmp[] = {
	-1,0,1,
	-2,0,2,
	-1,0,1
}; Mat<> SobelKernel(3, 3, SobelKernelTmp);
/*************************************************************************************************
*								边缘检测
*	[目的]: 标识数字图像中亮度变化明显的点.
*	[公式]: EdgeImage = Conv(Image , SobelKernel)
*************************************************************************************************/
Mat<>& EdgeDetection(Mat<>& in, Mat<>& out) {
	Mat<> out_x, out_y;
	out_x.conv(in, SobelKernel,                  1);
	out_y.conv(in, SobelKernel.transpose(out_y), 1);
	out.zero(in);
	for (int i = 0; i < in.size(); i++)
		out[i] = sqrt(out_x[i] * out_x[i] + out_y[i] * out_y[i]);
	return out;
}
/*************************************************************************************************
*								傅里叶变换
*	[目的]: 转频域图像.
*************************************************************************************************/
Mat<>& FourierTransform(Mat<>& in, Mat<>& out) {
	return out;
}
Mat<>& InvFourierTransform(Mat<>& in, Mat<>& out) {
	return out;
}
/*************************************************************************************************
*								Gauss 滤波
* [输入]: in: 输入原图 out: 输出图像  size: 核的大小  sigma: 正态分布标准差
*************************************************************************************************/
Mat<>& GaussFilter(Mat<>& in, int size, float sigma, Mat<>& out) {
	if (size <= 0 || sigma == 0) return out;
	Mat<> GaussKernel(size, size);
	for (int y = 0; y < size; y++)
		for (int x = 0; x < size; x++)
			GaussKernel(x, y) = 1 / (2 * PI * sigma * sigma) * exp(-(pow(x - size / 2, 2) + pow(y - size / 2, 2)) / (2 * sigma * sigma));
	return out.conv(in, GaussKernel *= 1 / GaussKernel.sum(), 1);
}
/*************************************************************************************************
*								直方图
* [目的]: 统计[0,255]亮度的像素个数分布.
*************************************************************************************************/
Mat<int>& Histograms(Mat<>& in, Mat<>& out) {
	out.zero(0xFF);
	for (int i = 0; i < in.size(); i++) out[(unsigned char)(in[i] * 0xFF)]++;
}

}
#endif