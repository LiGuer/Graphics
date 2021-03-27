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
#include "../LiGu_AlgorithmLib/Mat.h"
#include "../LiGu_AlgorithmLib/BasicMachineLearning.h"
#include "../opencv2-include/opencv2/opencv.hpp"
#pragma comment(lib,"../opencv2-include/opencv_world430.lib")

namespace DigitalImageProcessing {
#define PI 3.141592653589
	/*--------------------------------[ ͼ������/��� ]--------------------------------*/
	Mat<double>* Input(const char* inputImgUrl, Mat<double>* data) {
		cv::Mat input = cv::imread(inputImgUrl, cv::IMREAD_COLOR);
		for (int k = 0; k < 3; k++) data[k].zero(input.rows, input.cols);

		for (int i = 0; i < input.rows * input.cols; i++) {
			cv::MatIterator_<cv::Vec3b> it = input.begin<cv::Vec3b>() + i;
			for (int k = 0; k < 3; k++)  data[2 - k][i] = (double)(*it)[k] / 255;
		}
		return data;
	}
	void Output(const char* outputImgUrl, Mat<double>* data) {
		unsigned char* output = (unsigned char*)calloc(data[0].cols * data[0].rows * 3, sizeof(unsigned char));
		for (int i = 0; i < data[0].rows * data[0].cols; i++)
			for (int k = 0; k < 3; k++)
				output[i * 3 + k] = (data[k])[i] * 255;
		FILE* fp = fopen(outputImgUrl, "wb");
		fprintf(fp, "P6\n%d %d\n255\n", data[0].cols, data[0].rows);	// дͼƬ��ʽ����ߡ��������ֵ
		fwrite(output, 1, data[0].cols * data[0].rows * 3, fp);	// дRGB����
		fclose(fp);
		free(output);
	}
	void Output(const char* outputImgUrl, Mat<double>& data) {
		unsigned char* output = (unsigned char*)calloc(data.cols * data.rows, sizeof(unsigned char));
		for (int i = 0; i < data.rows * data.cols; i++)
				for (int k = 0; k < 3; k++)
					output[i] = data[i] * 255;
		FILE* fp = fopen(outputImgUrl, "wb");
		fprintf(fp, "P5\n%d %d\n255\n", data.cols, data.rows);	// дͼƬ��ʽ����ߡ��������ֵ
		fwrite(output, 1, data.cols * data.rows, fp);	// дRGB����
		fclose(fp);
		free(output);
	}
	/*--------------------------------[ ��ֵ�� ]--------------------------------
	[Ŀ��]: ������ֵ����ͼ���Ϊ���ڰ�ͼ.
	[��ʽ]: BinarizationImage = Image > threshold ? 1 : 0;
	---------------------------------------------------------------------------*/
	Mat<double>& Binarization(Mat<double>& input, Mat<double>& output, double threshold = 0.5) {
		output.zero(input.rows, input.cols);
		for (int i = 0; i < input.rows * input.cols; i++)
			output[i] = input[i] > threshold ? 1 : 0;
		return output;
	}
	/*--------------------------------[ ��ɫ���� ]--------------------------------
	[Ŀ��]: �򻯾���ͼ���е�ɫ��.
	[�㷨]: K-Mean��ֵ����
	-----------------------------------------------------------------------------*/
	Mat<double>* ColorCluster(Mat<double>* input, Mat<double>* output, int K = 3, int TimesMax = 0x7FFFFFFF) {
		// Process input & output
		Mat<double> data(3, input[0].rows * input[0].cols);
		for (int k = 0; k < 3; k++)
			for (int i = 0; i < input[0].rows; i++)
				for (int j = 0; j < input[0].cols; j++)
					data(k, i * input[0].cols + j) = (input[k])(i, j);
		for (int k = 0; k < 3; k++) output[k].zero(input[0].rows, input[0].cols);
		// Color Cluster
		time_t now; srand((unsigned)time(&now));
		Mat<double> Center;
		Mat<int> Cluster, Cluster_Cur;
		K_Mean(data, K, Center, Cluster, Cluster_Cur, TimesMax);
		for (int i = 0; i < K; i++)
			for (int j = 0; j < Cluster_Cur[i]; j++)
				for (int dim = 0; dim < 3; dim++)
					(output[dim])(Cluster(i, j) / output[0].cols, Cluster(i, j) % output[0].cols) = Center(dim, i);
		return output;
	}
	/*--------------------------------[ ��Ե��� ]--------------------------------
	[Ŀ��]: ��ʶ����ͼ�������ȱ仯���Եĵ�.
	[��ʽ]: EdgeImage = Conv(Image , SobelKernel)
	-----------------------------------------------------------------------------*/
	Mat<double>& EdgeDetection(Mat<double>& input, Mat<double>& output) {
		Mat<double> SobelKernel(3, 3);
		{
			double t[] = {
				-1,0,1,
				-2,0,2,
				-1,0,1
			}; SobelKernel.getData(t);
		}
		Mat<double> output_x, output_y;
		output_x.conv(input, SobelKernel, 1);
		output_y.conv(input, SobelKernel.transposi(output_y), 1);
		output.zero(input.rows, input.cols);
		for (int i = 0; i < input.rows * input.cols; i++)
			output[i] = sqrt(output_x[i] * output_x[i] + output_y[i] * output_y[i]);
		return output;
	}
	/*--------------------------------[ ����Ҷ�任 ]--------------------------------
	[Ŀ��]: תƵ��ͼ��.
	-------------------------------------------------------------------------------*/
	Mat<double>& FourierTransform(Mat<double>& input, Mat<double>& output) {
		return output;
	}
	Mat<double>& InvFourierTransform(Mat<double>& input, Mat<double>& output) {
		return output;
	}
	/*--------------------------------[ Gauss �˲� ]--------------------------------
	* [����]: input: ����ԭͼ dst: ģ��ͼ��  size: �˵Ĵ�С  sigma: ��̬�ֲ���׼��
	--------------------------------------------------------------------------------*/
	Mat<double>& GaussFilter(Mat<double>& input, int size, float sigma, Mat<double>& output) {
		if (size <= 0 || sigma == 0)return;
		//��άGauss������
		Mat<double> kernel(size, size);
		double sum = 0;
		for (int y = 0; y < size; y++) {
			for (int x = 0; x < size; x++) {
				kernel(x, y) = (1 / (2 * PI * sigma * sigma))
					* exp(-((x - size / 2) * (x - size / 2) + (y - size / 2) * (y - size / 2)) / (2 * sigma * sigma));
				sum += kernel(x, y);
			}
		}
		kernel.mult(1 / sum, kernel);
		//Gauss���
		output.conv(input, kernel, 1);
		return output;
	}
	/*--------------------------------[ ת�Ҷ�ͼ ]--------------------------------
	[Ŀ��]: RGB��ͨ���ϲ�Ϊ�Ҷ�һͨ��
	[��ʽ]: Gray = 0.3 R + 0.59 G + 0.11 B
	-----------------------------------------------------------------------------*/
	Mat<double>& Gray(Mat<double>* input, Mat<double>& output, double Rk = 0.3, double Gk = 0.59, double Bk = 0.11) {
		output.zero(input[0].rows, input[0].cols);
		Mat<double> tmp;
		output.add(output, tmp.mult(Rk / (Rk + Gk + Bk), input[0]));
		output.add(output, tmp.mult(Gk / (Rk + Gk + Bk), input[1]));
		output.add(output, tmp.mult(Bk / (Rk + Gk + Bk), input[2]));
		return output;
	}
	/*--------------------------------[ ֱ��ͼ ]--------------------------------
	* [Ŀ��]: ͳ��[0,255]���ȵ����ظ����ֲ�.
	---------------------------------------------------------------------------*/
	Mat<int>& Histograms(Mat<double>& input, Mat<double>& output) {
		output.zero(255, 1);
		for (int i = 0; i < input.rows; i++)
			for (int j = 0; j < input.cols; j++)
				output[(unsigned char)(input(i, j) * 255)]++;
	}
	Mat<int>& Histograms(Mat<double>* input, Mat<double>& output) {
		output.zero(255, 3);
		for (int k = 0; k < 3; k++)
			for (int i = 0; i < input->rows; i++)
				for (int j = 0; j < input->cols; j++)
					output((unsigned char)(input[k](i, j) * 255), k)++;
	}
	/*--------------------------------[ ���� ]--------------------------------
	[Ŀ��]: ������ɫ�����䲹ɫ
	[��ʽ]: InvImage = 1 - Image
	-------------------------------------------------------------------------*/
	Mat<double>& Invert(Mat<double>& input, Mat<double>& output) {
		output.mult(1, input.negative(output)); return output;
	}
	Mat<double>* Invert(Mat<double>* input, Mat<double>* output) {
		for (int k = 0; k < 3; k++)  Invert(input[k], output[k]); return output;
	}
}
/*//Example
int main() {
	Mat<double> img[3];
	DigitalImageProcessing::Input("IMG01.jpg", img);
	Mat<double> colorclusterImg[3];
	DigitalImageProcessing::ColorCluster(img, colorclusterImg, 3, 20);
	Mat<double> grayImg;
	DigitalImageProcessing::Gray(colorclusterImg, grayImg);
	Mat<double> edgedetectionImg;
	DigitalImageProcessing::EdgeDetection(grayImg, edgedetectionImg);
	DigitalImageProcessing::Output("img.ppm", edgedetectionImg);
}
*/
#endif