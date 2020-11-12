#ifndef _MAT_H
#define _MAT_H
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
template<class T>
class Mat
{
public:
/******************************************************************************
*                    基础数据
******************************************************************************/
	T* data = NULL;
	int rows = 0, cols = 0;
	/*---------------- 构造析构函数 ----------------*/
	Mat() { ; }
	Mat(const int _rows, const int _cols) { zero(_rows, _cols); }
	Mat(const Mat& a) { assign(a); }
	~Mat() { free(data); }
	/*---------------- 基础函数 ----------------*/
	void clean() {memset(data, 0, sizeof(T) * rows * cols);}		//清零 
	void error() { exit(-1);}
/******************************************************************************
*                    基础矩阵
******************************************************************************/
	/*---------------- 零元 ----------------*/
	void zero(const int _rows, const int _cols) {
		if (data != NULL)free(data);
		data = (T*)malloc(sizeof(T) * _rows * _cols);
		memset(data, 0, sizeof(T) * _rows * _cols);
		rows = _rows;	cols = _cols;
	}
	/*---------------- 单位元 ----------------*/
	void E(const int _rows) {
		zero(_rows, _rows);
		for (int i = 0; i < rows; i++) {
			data[i * cols + i] = 1;
		}
	}
/******************************************************************************
*                    运算
*	* "[]"取元素
	* max/min
	* [ = ]assign	* [ + ]add		* [ * ]mult 矩阵乘/数乘
******************************************************************************/
	/*---------------- "[]"取元素 ----------------*/
	T& operator[](int i) { return data[i]; }
	/*---------------- max/min ----------------*/
	T max() const {
		T maxdata = *data;
		for (int i = 0; i < rows * cols; i++)maxdata = maxdata >= data[i] ? maxdata : data[i];
		return maxdata;
	}
	T max(int index) {
		T maxdata = *data;
		for (int i = 0; i < rows * cols; i++)
			if (maxdata > data[i]) { maxdata = data[i]; index = i; }
		return maxdata;
	}
	T min() const {
		T mindata = *data;
		for (int i = 0; i < rows * cols; i++)mindata = mindata <= data[i] ? mindata : data[i];
		return mindata;
	}
	/*----------------赋矩阵 [ = ]----------------*/
	void assign(const Mat& a) {
		if (a.data == NULL)error();
		zero(a.rows, a.cols);
		memcpy(data, a.data, sizeof(T) * a.rows * a.cols);
	}
	/*----------------加法 [ + ]----------------*/
	void add(Mat& a, Mat& b, Mat& ans) {
		if (a.rows != b.rows || a.cols != b.cols)error();
		Mat ansTemp(a);
		for (int i = 0; i < a.rows * a.cols; i++)ansTemp[i] += b[i];
		// Save Ans
		if (ans.data != NULL)free(ans.data); ans.data = ansTemp.data; ansTemp.data = NULL;
		ans.rows = ansTemp.rows; ans.cols = ansTemp.cols;
	}
	/*----------------乘法 [ * ]----------------*/
	void mult(const Mat& a, const Mat& b, Mat& ans) {
		if (a.cols != b.rows) error();
		Mat ansTemp(a.rows, b.cols);
		for (int i = 0; i < a.rows; i++) {
			for (int j = 0; j < b.cols; j++) {
				T sum;
				memset(&sum, 0, sizeof(sum));
				for (int k = 0; k < a.cols; k++) {
					T aV = a.data[i * a.cols + k];
					T bV = b.data[k * b.cols + j];
					sum += aV * bV;
				}
				ansTemp.data[i * ansTemp.cols + j] = sum;
			}
		}
		// Save Ans
		if (ans.data != NULL)free(ans.data); ans.data = ansTemp.data; ansTemp.data = NULL;
		ans.rows = ansTemp.rows; ans.cols = ansTemp.cols;
	}
	void mult(const double a, const Mat& b, Mat& ans) {
		Mat ansTemp(b.rows, b.cols);
		for (int i = 0; i < b.rows * b.cols; i++)
			ansTemp.data[i] = a * b.data[i];
		// Save Ans
		if (ans.data != NULL)free(ans.data); ans.data = ansTemp.data; ansTemp.data = NULL;
		ans.rows = ansTemp.rows; ans.cols = ansTemp.cols;
	}
	/*----------------元素求和 [ sum() ]----------------*/
	void sum(int dim, Mat& ans) {
		int _col = 1, _row = 1;
		if (dim == 0)_row = rows;
		else if (dim == 1)_col = cols;
		ans.zero(_row, _col);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				ans[i] += data[i * cols + j];
			}
		}
	}
	/*----------------转置 [ trans() ]----------------*/
	void transposi(Mat& ans) {
		Mat ansTemp(cols, rows);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				ansTemp.data[j * rows + i] = data[i * cols + j];
			}
		}
		// Save Ans
		if (ans.data != NULL)free(ans.data); ans.data = ansTemp.data; ansTemp.data = NULL;
		ans.rows = ansTemp.rows; ans.cols = ansTemp.cols;
	}
	/*----------------余子式 [ comi ]----------------*/
	/*----------------取逆 [ inv ]----------------*/
	/*----------------行列式 [ abs() ]----------------*/
	/*--------------伴随矩阵 [ adj() ]----------------*/
	/*----------------特征值特征向量 [ eig() ]----------------
	*	特征方程: AX = λX
	*		A: 目标矩阵		X: 特征向量		λ: 特征值
	*	性质:
	*		若 R 为正交矩阵 (R'R = E),有B = R`¹A R , 使得 BY = λY, 特征值不变.
	*				又有 X = R Y.
	*	[算法]雅可比迭代:
	*	* 原理:
	*		对于目标实矩阵A, 构造正交矩阵序列 R1, R2, ... , Rn，
	*			D0 = A
	*			Dj = RjT Dj-1 Rj
	*			=> limj->∞ Dj = D = diag(λ1, λ2, ... , λn)
	*		当非对角元素接近0时，算法即可停止。
		*
			AR, 右乘只改变 pth col and qth col
				djp = c ajp - s ajq
				djq = s ajp + c ajq
			RA, 左乘只改变 pth row and qth row
			R'AR:
				djp = c ajp - s ajq
				djq = s ajp + c ajq
				dpp = c² app + s² aqq - 2 c s apq
				dqq = s² app + c² aqq + 2 c s apq
				dpq = ( c² - s² ) apq + c s ( app - aqq )
				其他元素对称性可得
		*	每一步使得非对角线 dpq dqp 为零
			对dpq: (c² - s²)/(cs) = (aqq - app)/apq
			令 s = sinΦ	c = cosΦ	t = tanΦ = s / c
			θ = cot(2Φ) = (aqq - qpp) / (2 * apq)
			tan2Φ = (aqq - qpp) / apq = 2 * tanΦ / (1 - tan²Φ)
			t² + 2tθ - 1 = 0
	*------------------------------------------------*/
	void eig(T esp, Mat& eigvec, Mat& eigvalue) {
		if (rows != cols)return;
		// init
		eigvalue.assign(*this);
		eigvec.E(rows);
		int n = rows;
		Mat<double> R, RT;
		// begin iteration
		while (true) {
			// Calculate row p and col q
			int p, q;
			T maxelement = eigvalue[1];
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					if (i != j && fabs(eigvalue[i * n + j]) >= maxelement) {
						maxelement = fabs(eigvalue[i * n + j]); p = i; q = j;
					}
				}
			}if (maxelement < esp)return;
			// eigvalue eigvec
			T theta = 0.5 * atan2(2 * eigvalue[p * n + q], eigvalue[q * n + q] - eigvalue[p * n + p]);
			T c = cos(theta), s = sin(theta);		// c,s
			R.E(n);
			R[p * n + p] = c; R[p * n + q] = s;	// R
			R[q * n + p] = -s; R[q * n + q] = c;
			R.transposi(RT);
			eigvalue.mult(RT, eigvalue, eigvalue);		// Dj = RjT Dj-1 Rj
			eigvalue.mult(eigvalue, R, eigvalue);
			eigvec.mult(eigvec, R, eigvec);				// X = R Y
		}
	}
};
#endif