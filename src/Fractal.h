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
#ifndef FRACTAL
#define FRACTAL
#include <math.h>
#include <complex>
#include <vector>
#include "../../LiGu_AlgorithmLib/Mat.h"
#include "GraphicsND.h"
#define PI 3.141692653589

namespace Fractal {
/*################################################################################################
isMandelbrotSet		(C, Z0, iterateTimes)
Mandelbrot			(C, min, max, resSize, iterateTimes, Set)
Julia				(C, min, max, resSize, iterateTimes, Set)
HilbertCurve_xy2d	(n, x, y)
HilbertCurve_d2xy	(n, distance, int& x, int& y)
PerlinNoise			(x, y, randomGridGradient)
PerlinNoise			(output, frequency)
FractalTree3D		(linesSt, linesEd, level, alpha, fork)
Boids				(birds)
BoidsRule			(birds, index)
################################################################################################*/
/*************************************************************************************************
*								Mandelbrot集
*	[定义]: Zn+1 = Zn² + C
			所有能使Zn+1不发散的复数点C, 所构成的集合,即 Mandelbrot Set
			(不发散,不一定收敛,有可能在几个不同点来回跳)
*	[性质]: |Zn|>2不可能收敛, 即Mandelbrot Set在半径为2的圆内.
*************************************************************************************************/
// 是否属于MandelbrotSet/JuliaSet
int isMandelbrotSet(std::complex<double> C, std::complex<double> Z0, int iterateTimes) {
	std::complex<double> z = Z0;
	for (int epoch = 0; epoch < iterateTimes; epoch++) {	// 迭代
		if (abs(z) > 2) return epoch;						// |Zn|>2不可能收敛
		z = z * z + C;										// Zn+1 = Zn² + C   
	} return 0;							//属于,输出0; 不属于,输出判断出的当次迭代数
}
void Mandelbrot(std::complex<double> min, std::complex<double> max, int resSize, int iterateTimes, Mat<int> Set) {
	Set.zero(resSize, resSize);
	double deltaReal = (max.real() - min.real()) / resSize,
	       deltaImag = (max.imag() - min.imag()) / resSize;
	for (int i = 0; i < resSize; i++) {
		for (int j = 0; j < resSize; j++) {
			double real = min.real() + deltaReal * i,
				   imag = min.imag() + deltaImag * j;
			Set(i, j) = isMandelbrotSet(
				std::complex<double>(real, imag), 
				std::complex<double>(0, 0), 
				iterateTimes
			);
		}
	}
}
/*************************************************************************************************
*								Julia集
*	[定义]: Zn+1 = Zn² + C
			对于某复数值C,所有能使Zn+1不发散的Z0的集合,即 Julia Set
			类似于. Mandelbrot 曼德布洛特集
*************************************************************************************************/
void Julia(std::complex<double> C, std::complex<double> min, std::complex<double> max, int resSize, int iterateTimes, Mat<int> Set) {
	Set.zero(resSize, resSize);
	double deltaReal = (max.real() - min.real()) / resSize,
	       deltaImag = (max.imag() - min.imag()) / resSize;
	for (int i = 0; i < resSize; i++) {
		for (int j = 0; j < resSize; j++) {
			double real = min.real() + deltaReal * i, 
				   imag = min.imag() + deltaImag * j;
			Set(i, j) = isMandelbrotSet(C, std::complex<double>(real, imag), iterateTimes);
		}
	}
}
/*************************************************************************************************
*								Hilbert 曲线
*	[定义]: 一种自相似的分形曲线
*	[生成方法]: 四象限复制四分，翻转左下、右下, 使左下末同左上初、右下初同右上末,能够最短连接.
*	[性质]:
		* 边长: nth 2^n
		* 长度: nth 2^n - 1 / 2^n
		* 因为是四进制自相似, 所以曲线上位置 distance, 可以不断判断子象限，按二进制在位上叠加
*	[用途]: 1D to 2D 的映射算法, 随 n 增加, 映射的位置趋于收敛
*	[1,2阶曲线]: "┌┐"	   "┌┐┌┐"    "┌┐┌┐ "
					︱︱  =>	︱︱︱︱ =>   ︱└┘︱
							┌┐┌┐      └┐┌┘
							︱︱︱︱       -┘└-
*	[程序]: rotation : 翻转坐标, 使Hilbert曲线性质, 自相似地适用于左下、右下象限
		[1] xy2d(): 坐标 -> 曲线位置    [2] d2xy(): 曲线位置 -> 坐标
*	[Author]: 1891.Hilbert
*	[Reference]: [1] Wikipedia.Hilbert curve
*************************************************************************************************/
int HilbertCurve_xy2d(int n, int x, int y) {	//convert (x,y) to d
	int distance = 0;
	for (int s = n / 2; s > 0; s /= 2) {
		int rx = (x & s) > 0, 
			ry = (y & s) > 0;
		distance += s * s * ((0b11 * rx) ^ ry);
		// rotation
		if (ry == 0) {
			if (rx == 1) { 
				x = n - 1 - x; 
				y = n - 1 - y;
			}
			int t = x; x = y; y = t;
		}
	}return distance;
}
void HilbertCurve_d2xy(int n, int distance, int& x, int& y) {//convert d to (x,y)
	x = y = 0;
	for (int s = 1; s < n; s *= 2) {
		int rx = 1 & (distance / 2), 
			ry = 1 & (distance ^ rx);
		// rotation
		if (ry == 0) {
			if (rx == 1) { 
				x = n - 1 - x; 
				y = n - 1 - y;
			}
			int t = x; x = y; y = t;
		}
		x += s * rx; 
		y += s * ry;
		distance /= 4;
	}
}
/*************************************************************************************************
*								Perlin Noise
*	Function to linearly interpolate between a0 and a1 , Weight w should be in the range [0.0, 1.0]
*	[流程]:
		[1] 格点随机梯度矢量
		[2] (x,y)与格点距离,梯度点积
		[3] 插值
*************************************************************************************************/
double PerlinNoise(double x, double y, Mat<>& randomGridGradient) {
	// 对四个格点
	int x0[] = { x, x + 1, x, x + 1 }, 
		y0[] = { y, y, y + 1, y + 1 };
	double n[4];
	for (int i = 0; i < 4; i++) {
		//[1][2]
		double random = randomGridGradient(x0[i], y0[i]);
		n[i] = (x - x0[i]) * cos(random) 
			 + (y - y0[i]) * sin(random);
	}
	//[3]
	double 
		sx = x - (int)x,
		sy = y - (int)y,
		ix0 = (n[1] - n[0]) * (3.0 - sx * 2.0) * sx * sx + n[0],
		ix1 = (n[3] - n[2]) * (3.0 - sx * 2.0) * sx * sx + n[2];
	return (ix1 - ix0) * (3.0 - sy * 2.0) * sy * sy + ix0;
}
Mat<>& PerlinNoise(Mat<>& output, int frequency) {
	Mat<> randomGridGradient;
	randomGridGradient.rands(frequency + 1, frequency + 1, 0, 256);
	for (int y = 0; y < output.cols; y++)
		for (int x = 0; x < output.rows; x++)
			output(x, y) = PerlinNoise(
				(double)x / output.rows * frequency,
				(double)y / output.cols * frequency,
				randomGridGradient
			);
	return output;
}
/*************************************************************************************************
*								三维分形树 Fractal Tree 3D
*************************************************************************************************/
void FractalTree3D(std::vector<Mat<>>& linesSt, std::vector<Mat<>>& linesEd, int level, double alpha, int fork = 3) {
	if (level <= 0) return;
	// 确定旋转矩阵
	Mat<> st = linesSt.back(), ed = linesEd.back(), direction, rotateAxis, rotateMat, zAxis(3), tmp; zAxis.set(0, 0, 1);
	direction.sub(ed, st);
	if (direction[0] != 0 || direction[1] != 0) {
		GraphicsND::rotate(
			rotateAxis.crossProduct(direction, zAxis),
			-acos(tmp.dot(direction, zAxis) / direction.norm()),
			tmp.zero(3), rotateMat.E(4)
		); rotateMat.block(1, 3, 1, 3, rotateMat);
	}
	else rotateMat.E(3);
	//递归
	double Lenth = direction.norm(); 
	Mat<> endPoint(3);
	for (int i = 0; i < fork; i++) {
		endPoint.set(
			sin(alpha) * cos((double)i * 2 * PI / fork), 
			sin(alpha) * sin((double)i * 2 * PI / fork), 
			cos(alpha)
		);
		endPoint.add(ed, endPoint.mul(0.7 * Lenth, endPoint.mul(rotateMat, endPoint)));
		linesSt.push_back(ed);
		linesEd.push_back(endPoint);
		FractalTree3D(linesSt, linesEd, level - 1, alpha);
	}
}
/******************************************************************************
*                    Boids 鸟群
*	[Rule]:
		[1] Collision Avoidance: avoid collisions with nearby flockmates
		[2] Velocity Matching  : attmpt to match velocity with nearby flockmates
		[3] Flock Centering    : attmpt to stay close to nearby flockmates
*	[过程]:
		[可见域内]:
			[1] 可见距离判定
			[2] 可见角度判定
		[Rule 1]: 可见域内，同各个个体的反向方向，除以距离，作为加速度 a1
		[Rule 2]: 可见域内, 质心速度, 作为加速度 a2
		[Rule 3]: 可见域内, 质心位置方向, 作为加速度 a3
		[障碍规避]:
			[1] 生成均匀球面点分布
			[2] 依次从速度方向,向四周,发出射线进行障碍检测,射线距离在可见域内
			[3] 选择第一束探测无障碍的射线, 作为避障加速度 a4
		[ v ]: v = (Σwi·ai)·dt
			   v 归一化单位矢量: ^v  =>  v = v0·^v
		[ r ]: r = v·dt = v0·dt·^v
*	[Referance]:
		[1] Thanks for https://github.com/SebLague/Boids
******************************************************************************/
struct BoidsBird { Mat<> r{ 3 }, v{ 3 }, a{ 3 }; };
//[ rule 规则 ]
void BoidsRule(std::vector<BoidsBird>& birds, int index, 
	double visualField, double visualAngle, double* weight, void(*obstacleAvoidance)(BoidsBird& birds)
) {
	Mat<> distance(3), avoidDirection(3), groupVelocity(3), groupCenter(3), tmp;
	int groupNum = 0;
	for (int i = 0; i < birds.size(); i++) {
		if ( i == index
		|| distance.add(birds[i].r,birds[index].r.negative(distance)).norm()                   > visualField
		|| distance.dot(distance,  birds[index].v) / (distance.norm() * birds[index].v.norm()) < visualAngle
		) continue;
		groupNum++;
		//[Rule 1]: collisionAvoid  [Rule 2]: velocityMatching  [Rule 3]: flockCentering
		avoidDirection+=(tmp.mul(1 / distance.norm(), distance));
		groupVelocity += birds[i].v;
		groupCenter   += distance;
	}
	avoidDirection *= -1;
	// Update Acceleration
	birds[index].a.zero(); if (groupNum == 0) return;
	birds[index].a += (avoidDirection.mul(weight[0], avoidDirection.normalized()));
	birds[index].a += (groupVelocity .mul(weight[1], groupVelocity .normalized()));
	birds[index].a += (groupCenter   .mul(weight[2], groupCenter   .normalized()));
	obstacleAvoidance(birds[index]);
}
//[ play 运行 ]
void Boids(std::vector<BoidsBird>& birds, void(*obstacleAvoidance)(BoidsBird& birds),
	double visualField, double visualAngle, double* weight, double dt = 1, double speed = 3 // 能见范围、能见角度cos、各规则权值
) {
	for (int i = 0; i < birds.size(); i++) BoidsRule(birds, i, visualField, visualAngle, weight, obstacleAvoidance);
	Mat<> tmp(3);
	for (int i = 0; i < birds.size(); i++) {
		birds[i].v += (tmp.mul(dt,         birds[i].a));
		birds[i].r += (tmp.mul(dt * speed, birds[i].v.normalized()));
	}
}
}
#endif