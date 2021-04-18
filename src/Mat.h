﻿/*
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
	T* data = NULL;													//数据堆叠方向: 行优先
	int rows = 0, cols = 0;
/******************************************************************************
*                    基础函数
-------------------------------------------------------------------------------
Mat();                                      //构造/析构函数
Mat(const int _rows, const int _cols);
Mat(const int _rows);
Mat(const Mat& a);
~Mat();
void error();                               //报错
int size();                                 //Size
Mat& fill(T a);                             //填充
void eatMat(Mat& a);                        //吃掉另一个矩阵(指针操作)
void swap(Mat& a);                          //交换数据 [ swap ]
******************************************************************************/
	/*---------------- 构造/析构函数 ----------------*/
	Mat() { ; }
	Mat(const int _rows, const int _cols) { zero(_rows, _cols); }
	Mat(const int _rows) { E(_rows); }
	Mat(const Mat& a) { *this = a; }
	~Mat() { delete data; }
	/*---------------- 报错  ----------------*/
	static void error() { exit(-1); }
	/*---------------- Size  ----------------*/
	inline int size() const { return rows * cols; }
	/*---------------- 填充  ----------------*/
	inline Mat& fill(T a) { for (int i = 0; i < size(); i++) data[i] = a; return *this; }
	/*---------------- 吃掉另一个矩阵(指针操作)  ----------------*/
	inline void eatMat(Mat& a) {
		if (data != NULL) delete data;
		data = a.data; a.data = NULL;
		rows = a.rows; cols = a.cols; a.rows = a.cols = 0;
	}
	/*----------------交换数据 [ swap ]----------------*/
	void swap(Mat& a) {
		T* tmp = a.data; a.data = data; data = tmp;
		int t = a.rows; a.rows = rows; rows = t;
		t = a.cols; a.cols = cols; cols = t;
	}
/******************************************************************************
*                    基础矩阵
-------------------------------------------------------------------------------
*	[0]alloc	分配空间   
	[1]zero		零元/清零   
	[2]E		单位元   
	[3]ones		全1元  
	[4]rands	随机元 
******************************************************************************/
	/*---------------- 分配空间 ----------------*/
	Mat& alloc(const int _rows, const int _cols) {
		if (_rows != rows || _cols != cols) {
			if (data != NULL) delete data;
			data = (T*)malloc(_rows * _cols * sizeof(T));
			rows = _rows; cols = _cols;
		} return *this;
	}
	/*---------------- 零元/清零 ----------------*/
	inline Mat& zero() { memset(data, 0, sizeof(T) * size()); return *this; }
	Mat& zero(const int _rows, const int _cols) {
		alloc(_rows, _cols); zero();  return *this;
	}
	/*---------------- 单位元 ----------------*/
	Mat& E(const int _rows) {
		zero(_rows, _rows);
		for (int i = 0; i < rows; i++) data[i * cols + i] = 1;
		return *this;
	}
	/*---------------- 全1元 ----------------*/
	Mat& ones(const int _rows, const int _cols) {
		alloc(_rows, _cols); fill(1); return *this;
	}
	/*---------------- 随机元 ----------------*/
	Mat& rands(const int _rows, const int _cols, T st, T ed) {
		alloc(_rows, _cols);
		for (int i = 0; i < size(); i++)
			data[i] = rand() / double(RAND_MAX) * (ed - st) + st;	//[st,ed)
		return *this;
	}
/******************************************************************************
*                    基础运算
-------------------------------------------------------------------------------
T& operator[](int i);                    // 索引元素
T& operator()(int i, int j);
T& operator()(int i);
T max();                                // max/min
T max(int& index);
T min();
T min(int& index);
bool operator==(const Mat& b);          //判断相等 [==]
Mat& operator=(const Mat& a);           //赋矩阵 [=] //不能赋值自己
Mat& getData(T* a);
Mat& getData(T x, T y);
Mat& getData(T x, T y, T z);
Mat& operator+=(Mat& a);                //加法 [add +]
Mat& add(Mat& a, Mat& b);
Mat& operator-=(Mat& a);                //减法 [sub -]
Mat& sub(Mat& a, Mat& b);
Mat& mult(Mat& a, Mat& b);              //乘法 [mult ×]
Mat& operator*=(const Mat& a);
Mat& operator*=(const double a);        //数乘 [mult ×]
Mat& mult(const double a, Mat& b);
Mat& divide(const double a, Mat& b);    //数除 [divide /]
static T dot(Mat& a, Mat& b);           //点乘 [dot ·]
T dot(Mat& a);
Mat& crossProduct(Mat& a, Mat& b);      //叉乘 [crossProduct ×]
Mat& crossProduct_(Mat& a, Mat& b);
Mat& elementMult(Mat& a, Mat& b);       //元素乘 [elementMult ×]
Mat& elementMult(Mat& a);
Mat& elementDivide(Mat& a, Mat& b);     //元素除 [elementDivide /]
Mat& elementDivide(Mat& a);
Mat& negative(Mat& ans);                //负 [negative -]
Mat& transposi(Mat& ans);               //转置 [transposi T]
T sum();                                //求和 [sum Σ]
static T sum(Mat& a);
Mat& sum(Mat& ans,int dim);
T product();                            //求积 [product Π]
T norm();                               //范数 [norm ||x||]
Mat& normalized();                      //归一化 [normalized]
T comi(int i0, int j0);                 //余子式 [comi]
Mat& inv(Mat& ans);                     //取逆 [inv x~¹]
T abs();                                //行列式 [abs |x|]
Mat& adjugate(Mat& ans);                //伴随矩阵 [adjugate A*]
void eig(T esp, Mat& eigvec, Mat& eigvalue);    //特征值特征向量 [eig]
Mat& solveEquations(Mat& b, Mat& x);            //解方程组 [solveEquations]
void LUPdecomposition(Mat& U, Mat& L, Mat<int>& P);                     //LUP分解 [LUPdecomposition]
Mat& diag(Mat& ans);                                                    //构造对角矩阵 [diag]
Mat<double>& conv(Mat& a, Mat& b, int padding = 0, int stride = 1);     //卷积 [conv]
-------------------------------------------------------------------------------
*	运算嵌套注意,Eg: b.add(b.mult(a, b), a.mult(-1, a));
		不管括号第一二项顺序,都是数乘,乘法,加法, 问题原因暂不了解，别用该形式。
* 	加减乘，即使是自己也不会影响，效率也不影响
******************************************************************************/
	/*---------------- "[]" "()"取元素 ----------------
	* 索引方向: 先纵再横.
	---------------------------------------------*/
	T& operator[](int i) { return data[i]; }
	T& operator()(int i, int j) { return data[i * cols + j]; }
	T& operator()(int i) { return data[i]; }
	/*---------------- max/min ----------------*/
	T max() const {
		T maxdata = *data;
		for (int i = 1; i < size(); i++)maxdata = maxdata >= data[i] ? maxdata : data[i];
		return maxdata;
	}
	T max(int& index) {
		T maxdata = *data; index = 0;
		for (int i = 1; i < size(); i++)
			if (maxdata < data[i]) { maxdata = data[i]; index = i; }
		return maxdata;
	}
	T min() const {
		T mindata = *data;
		for (int i = 1; i < size(); i++)mindata = mindata <= data[i] ? mindata : data[i];
		return mindata;
	}
	T min(int& index) {
		T mindata = *data; index = 0;
		for (int i = 1; i < size(); i++)
			if (mindata > data[i]) { mindata = data[i]; index = i; }
		return mindata;
	}
	/*----------------判断相等 [ == ]----------------*/
	bool operator==(const Mat& b) {
		if (rows != b.rows || cols != b.cols)return false;
		for (int i = 0; i < size(); i++)
			if (data[i] != b.data[i])return false;
		return true;
	}
	/*----------------赋矩阵 [ = ]----------------*/ //不能赋值自己
	Mat& operator=(const Mat& a) {
		if (a.data == NULL)error();
		alloc(a.rows, a.cols);
		memcpy(data, a.data, sizeof(T) * size());
		return *this;
	}
	Mat& getData(T* a) {
		memcpy(data, a, sizeof(T) * size());
		return *this;
	}
	Mat& getData(T x, T y) {
		if (rows != 2 || cols != 1)error();
		data[0] = x;
		data[1] = y;
		return *this;
	}
	Mat& getData(T x, T y, T z) {
		if (rows != 3 || cols != 1)error();
		data[0] = x;
		data[1] = y;
		data[2] = z;
		return *this;
	}
	/*----------------加法 [ add + ]----------------*/
	Mat& operator+=(Mat& a) {
		if (a.rows != rows || a.cols != cols)error();
		for (int i = 0; i < a.size(); i++)data[i] += a[i];
		return *this;
	}
	Mat& add(Mat& a, Mat& b) {
		if (a.rows != b.rows || a.cols != b.cols)error();
		alloc(a.rows, a.cols);
		for (int i = 0; i < a.size(); i++)data[i] = a[i] + b[i];
		return *this;
	}
	/*----------------减法 [ sub - ]----------------*/
	Mat& operator-=(Mat& a) {
		if (a.rows != rows || a.cols != cols)error();
		for (int i = 0; i < a.size(); i++)data[i] -= a[i];
		return *this;
	}
	Mat& sub(Mat& a, Mat& b) {
		if (a.rows != b.rows || a.cols != b.cols)error();
		alloc(a.rows, a.cols);
		for (int i = 0; i < a.size(); i++)data[i] = a[i] - b[i];
		return *this;
	}
	/*----------------乘法 [ mult × ]----------------*/
	Mat& mult(Mat& a, Mat& b) {
		if (a.cols != b.rows) error();
		Mat ansTmp(a.rows, b.cols);
		for (int i = 0; i < a.rows; i++)
			for (int j = 0; j < b.cols; j++) 
				for (int k = 0; k < a.cols; k++) 
					ansTmp(i, j) += a(i, k) * b(k, j);
		eatMat(ansTmp); return *this;
	}
	Mat& operator*=(const Mat& a) {
		if (cols != a.rows) error();
		Mat ansTmp(rows, a.cols);
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < a.cols; j++) 
				for (int k = 0; k < cols; k++)
					ansTmp(i, j) += (*this)(i, k) * a(k, j);
		eatMat(ansTmp); return *this;
	}
	/*----------------数乘 [ mult × ]----------------*/
	Mat& operator*=(const double a) {
		for (int i = 0; i < size(); i++) data[i] *= a; return *this;
	}
	Mat& mult(const double a, Mat& b) {
		alloc(b.rows, b.cols);
		for (int i = 0; i < size(); i++) data[i] = a * b[i];
		return *this;
	}
	/*----------------数除 [ divide / ]----------------*/
	Mat& divide(const double a, Mat& b) {
		alloc(b.rows, b.cols);
		for (int i = 0; i < size(); i++) data[i] = a / b[i];
		return *this;
	}
	/*----------------点乘 [ dot · ]----------------
	*	a·b = Σ ai·bi = aT * b
	**------------------------------------------------*/
	static T dot(Mat& a, Mat& b) {
		if (a.rows != b.rows || a.cols != b.cols) error();
		T ans = a[0] * b[0];
		for (int i = 1; i < a.size(); i++)ans += a[i] * b[i];
		return ans;
	}
	T dot(Mat& a) {
		if (a.rows != rows && a.cols != cols) error();
		T ans = data[0] * a[0];
		for (int i = 1; i < size(); i++)ans += data[i] * a[i];
		return ans;
	}
	/*----------------叉乘 [ crossProduct × ]----------------
	//####################### 暂时只三维
	*	𝑎 × 𝑏 ⃑ = | 𝑥		𝑦	 𝑧  |
					| 𝑥𝑎	𝑦𝑎	 za |
					| 𝑥𝑏	𝑦𝑏	 zb |
	**------------------------------------------------*/
	Mat& crossProduct(Mat& a, Mat& b) {
		if (a.rows != b.rows)error();
		Mat ansTmp(a.rows, a.cols);
		ansTmp[0] = a[1] * b[2] - a[2] * b[1];
		ansTmp[1] = a[2] * b[0] - a[0] * b[2];
		ansTmp[2] = a[0] * b[1] - a[1] * b[0];
		eatMat(ansTmp); return *this;
	}
	Mat& crossProduct_(Mat& a, Mat& b) {
		if (a.rows != b.rows)error();
		alloc(a.rows, a.cols);
		data[0] = a[1] * b[2] - a[2] * b[1];
		data[1] = a[2] * b[0] - a[0] * b[2];
		data[2] = a[0] * b[1] - a[1] * b[0];
		return *this;
	}
	/*----------------元素乘 [ elementMult × ]----------------
	**------------------------------------------------*/
	Mat& elementMult(Mat& a, Mat& b) {
		if (a.rows != b.rows || a.cols != b.cols) error();
		alloc(a.rows, a.cols);
		for (int i = 0; i < size(); i++)data[i] = a[i] * b[i];
		return*this;
	}
	Mat& elementMult(Mat& a) {
		if (rows != a.rows || cols != a.cols) error();
		for (int i = 0; i < size(); i++)data[i] *= a[i];
		return *this;
	}
	/*----------------元素除 [ elementDivide / ]----------------
	**------------------------------------------------*/
	Mat& elementDivide(Mat& a, Mat& b) {
		if (a.rows != b.rows || a.cols != b.cols) error();
		alloc(a.rows, a.cols);
		for (int i = 0; i < size(); i++)data[i] = a[i] / b[i];
		return*this;
	}
	Mat& elementDivide(Mat& a) {
		if (rows != a.rows || cols != a.cols) error();
		for (int i = 0; i < size(); i++)data[i] /= a[i];
		return *this;
	}
	/*----------------负 [ negative - ]----------------*/
	Mat& negative(Mat& ans) {
		ans.alloc(rows, cols);
		for (int i = 0; i < size(); i++) ans[i] = -data[i];
		return ans;
	}
	/*----------------转置 [ transposi T ]----------------*/
	Mat& transposi(Mat& ans) {
		Mat ansTmp(cols, rows);
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
				ansTmp(j, i) = data(i, j);
		ans.eatMat(ansTmp); return ans;
	}
	/*----------------求和 [ sum Σ ]----------------*/
	T sum() {
		T ans = data[0];
		for (int i = 1; i < size(); i++)ans += data[i];
		return ans;
	}
	static T sum(Mat& a) {
		T ans = a[0];
		for (int i = 1; i < size(); i++)ans += a[i];
		return ans;
	}
	Mat& sum(Mat& ans,int dim) {
		if (dim == 0) {				//对每一列求和
			Mat ansTmp(1, cols);
			for (int j = 0; j < cols; j++)
				for (int i = 0; i < rows; i++)
					ansTmp[j] += (*this)(i, j);
			ans.eatMat(ansTmp); return ans;
		}
		if (dim == 1) {				//对每一行求和
			Mat ansTmp(rows, 1);
			for (int i = 0; i < rows; i++)
				for (int j = 0; j < cols; j++)
					ansTmp[i] += (*this)(i, j);
			ans.eatMat(ansTmp); return ans;
		}
		error(); return ans;
	}
	/*----------------求积 [ product Π ]----------------*/
	T product() {
		T ans = data[0];
		for (int i = 1; i < size(); i++)ans *= data[i];
		return ans;
	}
	/*----------------范数 [ norm ||x|| ]----------------
	*	||a|| = sqrt(a·a)
	**-------------------------------------------*/
	T norm() { return sqrt(dot(*this, *this)); }
	/*----------------归一化 [ normalized ]----------------
	*	[定义]: 使得|| x || = 1
	**------------------------------------------------------*/
	Mat& normalized() {
		T t = norm();
		if (t == 0)return *this;
		for (int i = 0; i < size(); i++)data[i] /= t;
		return *this;
	}
	/*----------------余子式 [ comi ]----------------
	*	Mij: A 去掉第i行，第j列
	**-----------------------------------------------*/
	T comi(int i0, int j0) {
		Mat tmp(rows - 1, cols - 1);
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				if (i == i0 || j == j0)continue;
				tmp(i < i0 ? i : i - 1, j < j0 ? j : j - 1) = data[i * cols + j];
			}
		}
		return tmp.abs();
	}
	/*----------------取逆 [ inv x~¹ ]----------------
	*	[定义]: A A~¹ = E
	*	[方法]: 利用不断解线性方程组，对每一列求解.
	**------------------------------------------*/
	Mat& inv(Mat& ans) {
		if (rows != cols)error();
		Mat tmp(rows, cols);
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
			for (int i = 0; i < rows; i++)tmp(i, k) = x[i];
		}
		ans.eatMat(tmp); return ans;
	}
	/*----------------行列式 [ abs |x| ]----------------
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
	/*--------------伴随矩阵 [ adjugate A* ]----------------
	*	[定义]: 伴随矩阵A* 由(i,j)代数余子式Aij构成
				 [ A00  ... ]
			A* = | A01  Aij |
				 [ A02  ... ]
	*	[性质]: A* A = |A|
	**---------------------------------------------*/
	Mat& adjugate(Mat& ans) {
		Mat ansTmp(rows, cols);
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
				ansTmp(i, j) = ((i + j) % 2 == 0 ? 1 : -1) * comi(i, j);
		ans.eatMat(ansTmp); return ans;
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
		eigvalue = (*this);
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
	/*----------------构造对角矩阵 [ diag ]----------------*/
	Mat& diag(Mat& ans) {
		Mat ansTmp;
		if (rows == cols) {
			ansTmp.alloc(rows, 1);
			for (int i = 0; i < rows; i++)ansTmp[i] = data[i * rows + i];
		}
		else if (rows == 1 || cols == 1) {
			int n = rows > cols ? rows : cols;
			ansTmp.alloc(n, n);
			for (int i = 0; i < n; i++)ansTmp(i, i) = data[i];
		}
		else error();
		ans.eatMat(ansTmp); return ans;
	}
	/*----------------卷积 [ conv ]----------------*/
	Mat<double>& conv(Mat& a, Mat& b, int padding = 0, int stride = 1) {
		Mat ansTmp((a.rows - b.rows + 2 * padding) / stride + 1, (a.cols - b.cols + 2 * padding) / stride + 1);
		// for each element of output
		for (int y = 0; y < ansTmp.cols; y++) {
			for (int x = 0; x < ansTmp.rows; x++) {
				// for each element of b
				for (int ky = 0; ky < b.cols; ky++) {
					for (int kx = 0; kx < b.rows; kx++) {
						// get the corresponding element of a
						int xt = -padding + x * stride + kx, yt = -padding + y * stride + ky;
						ansTmp(x, y) += (xt < 0 || xt >= a.rows || yt < 0 || yt >= a.cols) ? 0 : a(xt, yt) * b(kx, ky);
					}
				}
			}
		}
		eatMat(ansTmp); return *this;
	}
/******************************************************************************
*                    特殊操作
-------------------------------------------------------------------------------
Mat& getCol(int _col, Mat& a)               //读/写一列 [getCol/setCol]
Mat& setCol(int _col, Mat& a)
Mat& getRow(int _row, Mat& a)               //读/写一行 [getRow/setRow]
Mat& block(int rowSt, int rowEd, int colSt, int colEd, Mat& ans)	//子矩阵 [block]
Mat& horizStack(Mat& a, Mat& b)             //水平向拼接 [horizStack ]
******************************************************************************/
	/*----------------读/写一列 [getCol/setCol]----------------*/
	Mat& getCol(int _col, Mat& a) {
		a.alloc(rows, 1);
		for (int i = 0; i < rows; i++) a[i] = *this(i, _col);
		return a;
	}
	Mat& getRow(int _row, Mat& a) {
		a.alloc(1, cols);
		for (int i = 0; i < cols; i++) a[i] = *this(_row, i);
		return a;
	}
	Mat& setCol(int _col, Mat& a) {
		for (int i = 0; i < rows; i++)data[i * cols + _col] = a[i];
		return a;
	}
	/*----------------子矩阵 [block]----------------*/
	Mat& block(int rowSt, int rowEd, int colSt, int colEd, Mat& ans) {
		Mat ansTmp(rowEd - rowSt + 1, colEd - colSt + 1);
		for (int i = rowSt; i <= rowEd; i++)
			for (int j = colSt; j <= colEd; j++)
				ansTmp(i - rowSt, j - colSt) = (*this)(i, j);
		ans.eatMat(ansTmp); return ans;
	}
	/*----------------水平向拼接 [horizStack]----------------*/
	Mat& horizStack(Mat& a, Mat& b) {
		if (a.rows != b.rows)error();
		Mat ansTmp(a.rows, a.cols + b.cols);
		for (int i = 0; i < ansTmp.row; i++)
			for (int j = 0; j < ansTmp.cols; j++)
				ansTmp.data[i * cols + j] = j < a.cols ? a(i, j) : b(i, j - a.cols);
		eatMat(ansTmp); return *this;
	}
	/*----------------复制拓展 [ repeatCol  ]----------------*/
	Mat& repeatCol(int repeatNum, Mat& ans) {
		if (cols != 1)error();
		Mat ansTmp(rows, repeatNum);
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < repeatNum; j++)
				ansTmp(i, j) = data[i];
		ans.eatMat(ansTmp); return ans;
	}
};
#endif