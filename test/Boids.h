#ifndef BOIDS_H
#define BOIDS_H
#include "LiGu_Graphics/Graphics3D.h"
#include <ctime> 
void Delay(int time){clock_t now = clock(); while (clock() - now < time); }//time*1000为秒数 
/******************************************************************************
*                    Boids
*	[Rule]:
		[1] Collision Avoidance: avoid collisions with nearby flockmates
		[2] Velocity Matching: attempt to match velocity with nearby flockmates
		[3] Flock Centering: attempt to stay close to nearby flockmates
*	[Referance]:
		[1] Thanks for https://github.com/SebLague/Boids
******************************************************************************/
class Boids {
public:
	struct Boids_cell { Mat<double> r{ 3,1 }, v{ 3,1 }, a{ 3,1 }; };
	Graphics3D g{ 1000, 1000 };
	double visualField = 30, visualAngle = -0.5;	// 能见范围 // 能见角度cosθ
	double deltaT = 1, speed = 3;
	double weight[4] = { 1,1,1,4 };					// 各规则的权值
	Mat<double> obstacleAvoidance_direction[300];
	/*********************************[ 函数 ]*********************************/
	/*--------------------------------[ play 运行 ]--------------------------------*/
	void play(Boids_cell cell[], int n) {
		for (int i = 0; i < n; i++) {
			rule(cell, n, i);
			obstacleAvoidance(cell[i]);
		}
		Mat<double> temp(3, 1);
		for (int i = 0; i < n; i++) {
			cell[i].v.add(cell[i].v, temp.mult(deltaT, cell[i].a));
			cell[i].v.normalization();
			cell[i].r.add(cell[i].r, temp.mult(deltaT * speed, cell[i].v));
		}
	}
	/*--------------------------------[ rule 规则 ]--------------------------------*/
	void rule(Boids_cell cell[], int n,int index) {
		Mat<double> distance(3, 1);
		Mat<double> avoidDirection(3, 1),groupVelocity(3, 1), groupCenter(3, 1),temp;
		int groupNum = 0;
		for (int i = 0; i < n; i++) {
			if (i == index)continue;
			// 能见范围
			distance.add(cell[i].r, cell[index].r.negative(distance));
			if (distance.norm() > visualField)continue;
			// 能见角度
			double angle = distance.dot(distance, cell[index].v) / (distance.norm() * cell[index].v.norm()); // cosθ = a * b / (|a| * |b|)
			if (angle < visualAngle)continue;
			groupNum++;
			//[Rule 1]: collisionAvoid
			avoidDirection.add(avoidDirection, temp.mult(1.0 / distance.norm(), distance.negative(temp)));
			//[Rule 2]: velocityMatching  [Rule 3]: flockCentering
			groupVelocity.add(groupVelocity, cell[i].v);
			groupCenter.add(groupCenter, distance);
		}
		// Update Acceleration
		cell[index].a.clean();
		if (groupNum == 0) return;
		avoidDirection.mult(weight[0] * 1.0 / groupNum, avoidDirection);
		groupVelocity.mult(weight[1] * 1.0 / groupNum, groupVelocity);
		groupCenter.mult(weight[2] * 1.0 / groupNum, groupCenter);
		cell[index].a.add(cell[index].a, avoidDirection.normalization());
		cell[index].a.add(cell[index].a, groupVelocity.normalization());
		cell[index].a.add(cell[index].a, groupCenter.normalization());
	}
	/*--------------------------------[ obstacleAvoidance 障碍规避 ]--------------------------------*/
	void obstacleAvoidance_init() {
		int numViewDirections = 300;
		// 均匀球面点
		double goldenRatio = (1 + sqrt(5)) / 2;				// 黄金分割点
		double angleIncrement = PI * 2 * goldenRatio;
		for (int i = 0; i < numViewDirections; i++) {
			double t = (double)i / numViewDirections, inclination = acos(1 - 2 * t), azimuth = angleIncrement * i;
			obstacleAvoidance_direction[i].zero(3, 1);
			obstacleAvoidance_direction[i][0] = sin(inclination) * cos(azimuth);
			obstacleAvoidance_direction[i][1] = sin(inclination) * sin(azimuth);
			obstacleAvoidance_direction[i][2] = cos(inclination);
		}
	}
	void obstacleAvoidance(Boids_cell& cell) {
		// 确定坐标变换矩阵
		Mat<double>rorateAxis(3, 1), rotateMat, Temp(3, 1);
		rorateAxis.crossProduct(cell.v, obstacleAvoidance_direction[0]);
		double theta = -acos(Temp.dot(cell.v, obstacleAvoidance_direction[0])); // cosθ = a * b / (|a| * |b|)
		g.rotate(rorateAxis.normalization(), theta, Temp, rotateMat);

		Temp = rotateMat; rotateMat.E(3);
		for (int i = 0; i < 3; i++)for (int j = 0; j < 3; j++)rotateMat(i, j) = Temp(i, j);
		// 碰撞检测
		Mat<double>derection(3, 1);
		for (int k = 0; k < 300; k++) {
			derection.mult(rotateMat, obstacleAvoidance_direction[k]);
			derection[0] = -derection[0];
			bool flag = 1; int distance = 1;
			while (distance < visualField) {
				Temp.add(Temp.mult(distance, derection), cell.r);
				if (obstacleFunction(Temp)) { flag = 0; break; }
				distance += 1;
			}
			if (flag)break;
		}
		cell.a.add(cell.a, derection.mult(weight[3], derection));
	}
	bool obstacleFunction(Mat<double>& point) {
		double size = 300;
		for (int i = 0; i < 3; i++) {
			if (fabs(point[i])>= 300)return true;
		}
		return false;
	}
	/*--------------------------------[ Boidss 主函数 ]--------------------------------*/
	void Boids_main() {
		double size = 300;
		Mat<double> Axis(3, 1); Axis[0] = 0.7; Axis[1] = 0.4; Axis[2] = 0.59;
		Mat<double> center(3, 1); center[0] = 500; center[1] = 500; center[2] = 500;
		double theta = 1;

		int N = 200;
		Boids_cell cell[200];
		for (int i = 0; i < N; i++) {
			cell[i].r.rands(3, 1, -size / 2, size / 2);
			cell[i].v.rands(3, 1, 0, 1);
			cell[i].v.normalization();
		}

		Mat<double> temp;
		Mat<double> zero0(3, 1), zero(3, 1), zero2(3, 1); zero2[0] = size; zero2[1] = size; zero2[2] = size; zero[0] = -size; zero[1] = -size; zero[2] = -size;
		obstacleAvoidance_init();
		while (1) {
			theta += 0.01; if (theta >= 2 * 3.1415926) theta = 0;
			g.TransformMat.E(4);
			g.rotate(Axis, theta, zero0);
			g.translation(center);

			play(cell, N);
			for (int i = 0; i < N; i++) {
				g.g->PaintColor = 0xFFFFFF;
				g.drawLine(cell[i].r, temp.add(temp.mult(10, cell[i].v), cell[i].r));
			}
			Delay(20);
			g.g->PicWrite("D:/LIGU.ppm");
			g.g->clear(0);
			for (int i = 0; i < 3; i++) {
				g.drawCuboid(zero, zero2);
			}
		}
	}
};

/*
int main()
{
	Delay(1000);
	Boids boids;
	boids.Boids_main();
}
*/
#endif
