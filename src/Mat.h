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

Reference.
[1]Introduction Algorithms.THOMAS H.CORMEN,CHARLES E.LEISERSON,RONALD L.RIVEST,CLIFFORD STEIN
==============================================================================*/
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
*                    核心数据
******************************************************************************/
	T* data = NULL;
	int rows = 0, cols = 0;
/******************************************************************************
*                    基础函数
******************************************************************************/
	/*---------------- 构造析构函数 ----------------*/
	Mat() { ; }
	Mat(const int _rows, const int _cols) { zero(_rows, _cols); }
	Mat(const Mat& a) { *this = a; }
	~Mat() { free(data); }
	/*---------------- 基础函数 ----------------*/
	void clean() {memset(data, 0, sizeof(T) * rows * cols);}		//清零 
	void error() { exit(-1);}
	void eatMat(Mat& a) {											//吃掉另一个矩阵的数据 (指针操作)
		if (data != NULL)free(data); 
		data = a.data; a.data = NULL;
		rows = a.rows; cols = a.cols; a.rows = a.cols = 0;
	}
/******************************************************************************
*                    基础矩阵
*	[1] 零元 zero		[2] 单位元 E		[3] 随机元 rands
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
	/*---------------- 随机元 ----------------*/
	void rands(const int _rows, const int _cols,T st,T ed) {
		zero(_rows, _cols);
		for (int i = 0; i < rows * cols; i++) {
			data[i] = rand() / double(RAND_MAX) * (ed - st) + st;	//[st,ed)
		}
	}
/******************************************************************************
*                    基础运算
-------------------------------------------------------------------------------
T& operator[](int i)                        // "[]"取元素
T& operator()(int i, int j)                 // "()"取元素
T& operator()(int i)
T max()                                     // max/min
T max(int& index)
T min()
T min(int& index)
Mat& operator=(const Mat& a)                //赋矩阵 [ = ]  //不能赋值自己
Mat& add(Mat& a, Mat& b)                    //加法 [ add ]
Mat& mult(const Mat& a, const Mat& b)       //乘法 [ mult ]
Mat& mult(const double a, const Mat& b)     //数乘 [ mult ]
Mat& dot(const Mat& a, const Mat& b)        //点乘 [ dot ]
Mat& negative(Mat& ans)                     //负 [ negative ]
Mat& transposi(Mat& ans)                    //转置 [ trans ]
void sum(int dim, Mat& ans)                 //元素求和 [ sum ]
T norm()                                    //范数 [ norm ]
T comi(int i0, int j0)                      //余子式 [ comi ]
Mat& inv(Mat& ans)                          //取逆 [ inv ]
T abs()                                     //行列式 [ abs ]
Mat& adjugate(Mat& ans)                     //伴随矩阵 [ adjugate ]
void eig(T esp, Mat& eigvec, Mat& eigvalue) //特征值特征向量 [ eig ]
Mat& solveEquations(Mat& b, Mat& x)         //解方程组 [ solveEquations ]
void LUPdecomposition(Mat& U, Mat& L, Mat& P) //LUP分解 [ LUPdecomposition ]
-------------------------------------------------------------------------------
*	运算嵌套注意,Eg: b.add(b.mult(a, b), a.mult(-1, a)); 
		不管括号第一二项顺序,都是数乘,乘法,加法, 问题原因暂不了解，别用该形式。
******************************************************************************/
	/*---------------- "[]"取元素 ----------------*/
	T& operator[](int i) { return data[i]; }
	T& operator()(int i, int j) { return data[i * cols + j]; }
	T& operator()(int i) { return data[i]; }
	/*---------------- max/min ----------------*/
	T max() const {
		T maxdata = *data;
		for (int i = 0; i < rows * cols; i++)maxdata = maxdata >= data[i] ? maxdata : data[i];
		return maxdata;
	}
	T max(int& index) {
		T maxdata = *data; index = 0;
		for (int i = 0; i < rows * cols; i++)
			if (maxdata < data[i]) { maxdata = data[i]; index = i; }
		return maxdata;
	}
	T min() const {
		T mindata = *data;
		for (int i = 0; i < rows * cols; i++)mindata = mindata <= data[i] ? mindata : data[i];
		return mindata;
	}
	T min(int& index) {
		T mindata = *data; index = 0;
		for (int i = 0; i < rows * cols; i++)
			if (mindata > data[i]) { mindata = data[i]; index = i; }
		return mindata;
	}
	/*----------------赋矩阵 [ = ]----------------*/ //不能赋值自己
	Mat& operator=(const Mat& a) {
		if (a.data == NULL)error();
		zero(a.rows, a.cols);
		memcpy(data, a.data, sizeof(T) * a.rows * a.cols);
		return *this;
	}
	/*----------------加法 [ add ]----------------*/
	Mat& add(Mat& a, Mat& b) {
		if (a.rows != b.rows || a.cols != b.cols)error();
		Mat ansTemp(a);
		for (int i = 0; i < a.rows * a.cols; i++)ansTemp[i] += b[i];
		eatMat(ansTemp);
		return *this;
	}
	/*----------------乘法 [ mult ]----------------*/
	Mat& mult(const Mat& a, const Mat& b) {
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
		eatMat(ansTemp);
		return *this;
	}
	/*----------------数乘 [ mult ]----------------*/
	Mat& mult(const double a, const Mat& b) {
		Mat ansTemp(b.rows, b.cols);
		for (int i = 0; i < b.rows * b.cols; i++)
			ansTemp.data[i] = a * b.data[i];
		eatMat(ansTemp);
		return *this;
	}
	/*----------------点乘 [ dot ]----------------
	*	a·b = Σ ai·bi = aT * b
	**------------------------------------------------*/
	T dot(const Mat& a, const Mat& b) {
		T ans;
		memset(ans, 0, sizeof(T));
		for (int i = 0; i < rows; i++)ans += a[i] * b[i];
		return ans;
	}
	/*----------------负 [ negative ]----------------*/
	Mat& negative(Mat& ans) {
		Mat ansTemp(*this);
		for (int i = 0; i < rows * cols; i++)
			ansTemp[i] = -ansTemp[i];
		ans.eatMat(ansTemp);
		return ans;
	}
	/*----------------转置 [ transposi ]----------------*/
	Mat& transposi(Mat& ans) {
		Mat ansTemp(cols, rows);
		for (int i = 0; i < rows; i++) 
			for (int j = 0; j < cols; j++) 
				ansTemp.data[j * rows + i] = data[i * cols + j];
		ans.eatMat(ansTemp);
		return ans;
	}
	/*----------------元素求和 [ sum ]----------------*///########
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
	/*----------------范数 [ norm ]----------------
	*	||a|| = sqrt(a·a)
	**-------------------------------------------*/
	T norm() { return sqrt(dot(*this, *this)); }
	/*----------------余子式 [ comi ]----------------
	*	Mij: A 去掉第i行，第j列
	**-----------------------------------------------*/
	T comi(int i0, int j0) {
		Mat temp(rows - 1, cols - 1);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				if (i == i0 || j == j0)continue;
				temp(i < i0 ? i : i - 1, j < j0 ? j : j - 1) = data[i * cols + j];
			}
		}
		return temp.abs();
	}
	/*----------------取逆 [ inv ]----------------
	*	[定义]: A A~¹ = E
	*	[方法]: 利用不断解线性方程组，对每一列求解.
	**------------------------------------------*/
	Mat& inv(Mat& ans) {
		if (rows != cols)error();
		Mat temp(rows, cols);
		int n = rows;
		// LUP分解
		Mat L, U; Mat<int> P;
		LUPdecomposition(U, L, P);
		//对每一列
		Mat b(n, 1), x(n, 1);
		for (int k = 0; k < n; k++) {
			b.clean(); b[k] = 1;
			// 解线性方程组
			//solve y
			for (int i = 0; i < n; i++) {
				x[i] = b[P[i]];		//yi
				for (int j = 0; j < i; j++) x[i] -= x[j] * L(i, j);
			}
			//solve x
			for (int i = n - 1; i >= 0; i--) {
				for (int j = i + 1; j < n; j++) x[i] -= x[j] * U(i, j);
				x[i] /= U(i, i);
			}
			//合并至结果
			for (int i = 0; i < rows; i++)temp(i, k) = x[i];
		}
		ans.eatMat(temp);
		return ans;
	}
	/*----------------行列式 [ abs ]----------------
	*	|A| = Σiorj aij·Aij
	*	Aij = (-1)^(i+j)·Mij		// Mij余子式
	**----------------------------------------------*/
	T abs() {
		if (rows != cols)error();
		if (rows == 1)return data[0];
		T ans;
		memset(&ans, 0, sizeof(T));
		for (int i = 0; i < rows; i++)
			ans += data[i * cols] * (i % 2 == 0 ? 1 : -1) * comi(i, 0);
		return ans;
	}
	/*--------------伴随矩阵 [ adjugate ]----------------
	*	[定义]: 伴随矩阵A* 由(i,j)代数余子式Aij构成
				 [ A00  ... ]
			A* = | A01  Aij |
			     [ A02  ... ]
	*	[性质]: A* A = |A|
	**---------------------------------------------*/
	Mat& adjugate(Mat& ans) {
		ans.zero(rows, cols);
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
				ans(i, j) = ((i + j) % 2 == 0 ? 1 : -1)* comi(i, j);
	}
	/*----------------特征值特征向量 [ eig ]----------------
	*	[定义]: 特征方程: AX = λX
	*		A: 目标矩阵		X: 特征向量		λ: 特征值
	*	[性质]:
	*		若 R 为正交矩阵 (R'R = E),有B = R~¹A R , 使得 BY = λY, 特征值不变.
	*				又有 X = R Y.
	*	[算法]: 雅可比迭代:
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
		//[1] init
		eigvalue.assign(*this);
		eigvec.E(rows);
		int n = rows;
		Mat<double> R, RT;
		//[2] begin iteration
		while (true) {
			//[3] Calculate row p and col q
			int p, q;
			T maxelement = eigvalue[1];
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					if (i != j && fabs(eigvalue[i * n + j]) >= maxelement) {
						maxelement = fabs(eigvalue[i * n + j]); p = i; q = j;
					}
				}
			}if (maxelement < esp)return;			// [2]
			//[4] eigvalue eigvec
			T theta = 0.5 * atan2(2 * eigvalue[p * n + q], eigvalue[q * n + q] - eigvalue[p * n + p]);
			T c = cos(theta), s = sin(theta);		// c,s
			R.E(n);
			R[p * n + p] = c; R[p * n + q] = s;		// R
			R[q * n + p] = -s; R[q * n + q] = c;
			R.transposi(RT);
			eigvalue.mult(RT, eigvalue);			// Dj = RjT Dj-1 Rj
			eigvalue.mult(eigvalue, R);
			eigvec.mult(eigvec, R);					// X = R Y
		}
	}
	/*----------------解方程组 [ solveEquations ]----------------
	*	[定义]: A x = b
	*			ps.直接x = b A~¹ 会存在数值不稳定现象
	*	[算法]: LUP分解
	*	[推导]
			P A = L U
			L: 单位下三角矩阵  U: 上三角矩阵  P: 置换矩阵
			*	置换矩阵:是一个方阵，每行和每列只有一个1，其他均为0，表示矩阵初等行变化
			*	置换矩阵是可逆的
			*	因为置换矩阵每行只有一个1，可以变为一维数组，每行计入改行1的位置
			*	P b 是对b的行互换，即.由上条，相当于b[Pi]
			可得 L U x = P b
			令 y = U x  =>  L y = P b	解得 y
			代入 U x = y  解得 x
			即. A x = P~¹ L U x = P~¹ L y = P~¹ P b = b
	*	[过程]:
			[1] LUP分解
			[2] LUP-Solve
				[3] Solve y:
					for i = 1 to n
								  j=1toi-1   
						yi = b[Pi] - Σ   Lij yj
				[4] Solve x:
					for i = n to 1
								j=i+1toN   
						x = ( yi - Σ  uij xj ) / uii
	**--------------------------------------------*/
	Mat& solveEquations(Mat& b, Mat& x) {
		int n = rows;
		x.zero(n, 1);
		//[1] LUP分解
		Mat U, L; Mat<int> P;
		LUPdecomposition(U, L, P);
		//[2] LUP - Solve
		//[3] solve y
		for (int i = 0; i < n; i++) {
			x[i] = b[P[i]];		//yi
			for (int j = 0; j < i; j++) x[i] -= x[j] * L(i, j);
		}
		//[4] solve x
		for (int i = n - 1; i >= 0; i--) {
			for (int j = i + 1; j < n; j++) x[i] -= x[j] * U(i, j);
			x[i] /= U(i, i);
		}
		return x;
	}
	/*----------------LUP分解 [ LUPdecomposition ]----------------
	*	[定义]: P A = L U		其中 L: 单位下三角矩阵  U: 上三角矩阵  P: 置换矩阵
			*	因为置换矩阵每行只有一个1，可以变为一维数组，每行计入改行1的位置
	*	[算法]: 高斯消元法
			[1] 从其他方程中减去第1方程的倍数，以把那些方程第1变量消去。
			[2] 从第3及以后方程中减去第2方程倍数，以把这些方程的第1,2变量都消去。
			[3] 重复过程，直至变为上三角矩阵U，单位下三角L是由消去变量所用行的乘数组成

			* 主元pivot: LPU分解中所除元素称为主元，它们处于矩阵U的对角线上。
			* 选主元: 采用置换避免除0，避免除数很小(数值会不稳定)的操作
			* 把第1行与第k行互换 <=> 置换矩阵Q左乘A--QA
	*	[过程]:
			[1] 对于每一列
				[2] 选主元
				[3] 置换行,记录在P中
				[4] LU分解: 高斯消元法
			[5] A中包含U,L，分离出来即可
	**---------------------------------------------*/
	void LUPdecomposition(Mat& U, Mat& L, Mat<int>& P) {
		if (rows != cols)error();
		int n = rows;
		Mat A(*this);
		P.zero(n, 1);
		for (int i = 0; i < n; i++)P[i] = i;
		//[1]
		for (int k = 0; k < n; k++) {
			//[2] 选主元
			T maxvalue = 0;
			int kt;
			for (int i = k; i < n; i++) {
				if (fabs(A(i, k)) > maxvalue) { maxvalue = fabs(A(i, k)); kt = i; }
			}
			if (maxvalue == 0)error();	// singular matrix，秩 rank<n
			//[3] 置换行
			for (int i = 0; i < n; i++) {
				T t = A(k, i); A(k, i) = A(kt, i); A(kt, i) = t;
			}
			int t = P[k]; P[k] = P[kt]; P[kt] = t;
			//[4] LU分解: 高斯消元法
			for (int i = k + 1; i < n; i++) {
				A(i, k) /= A(k, k);		//aik存储消去该行第k位所需的乘数,即L
				for (int j = k + 1; j < n; j++)
					A(i, j) -= A(i, k) * A(k, j);	//初等行变换，消去该行第k位
			}
		}
		//[5] A中包含U,L，分离出来即可
		U.zero(n, n); L.E(n);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i > j)L(i, j) = A(i, j);
				else U(i, j) = A(i, j);
			}
		}
	}
/******************************************************************************
*                    特殊操作
******************************************************************************/
	/*----------------水平向拼接 [ horizStack ]----------------*/
	Mat& horizStack(Mat& a, Mat& b) {
		if (a.rows != b.rows)error();
		Mat ansTemp(a.rows, a.cols + b.cols);
		for (int i = 0; i < ansTemp.row; i++)
			for (int j = 0; j < ansTemp.cols; j++)
				ansTemp.data[i * cols + j] = j < a.cols ? a(i, j) : b(i, j - a.cols);
		eatMat(ansTemp);
		return *this;
	}
	/*----------------交换数据 [ swap ]----------------*/
	void swap(Mat& a) {
		T* tptr = a.data;a.data = data;data = tptr;
		int t = a.rows; a.rows = rows; rows = t;
		t = a.cols; a.cols = cols; cols = t;
	}
	/*----------------得到一列 [ getCol ]----------------*/
	Mat& getCol(int _col, Mat& a) {
		a.zero(rows, 1);
		for (int i = 0; i < rows; i++) a[i] = data[i * cols + _col];
		return a;
	}
};
#endif