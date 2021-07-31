# LiGu-Graphics 计算机图形学  
* <Graphics.h/cpp>				图形学2D-Pix
* <GraphicsND.h/cpp>			任意维图形学
* <Plot.h/cpp>					科学制图
* <RayTracing.h/cpp>			光线追踪
* <RGB.h>						像素颜色
* <Fractal.h>					分形
* <ComputationalGeometry.h>		计算几何
* <DigitalImageProcessing.h>	数字图像处理
* <GraphicsFileCode.h>			图形文件编译码
* <ReadImg.exe>					实时动态显示图片

## API
### <Graphics.h> 图形学2D-Pix:  
```
/*-------------------------------- 基础参数 --------------------------------*/
Mat<RGB>	Canvas{ 100, 100 };											//图
Mat<FP64>	TransMat;													//变换矩阵
ARGB PaintColor = 0xFFFFFF;												//画笔颜色
INT32S 
	PaintSize = 0,														//画笔大小
	FontSize  = 16;														//字符大小
/*-------------------------------- 底层函数 --------------------------------*/
Graphics() { ; }
~Graphics() { }															//析构函数
Graphics (INT32S width, INT32S height) { init(width, height); }	
void init(INT32S width = 100, INT32S height = 100);						//初始化
void clear(ARGB color);	 												//清屏
void setPoint		(INT32S x, INT32S y, ARGB color);					//底层画点
ARGB readPoint		(INT32S x, INT32S y); 								//读点 
void readImg	(const char* filename);									//读图
void writeImg	(const char* filename);									//存图
bool judgeOutRange	(INT32S x0, INT32S y0);								//判断过界
void transSelf();														//全图变换
void CutSelf		(INT32S sx, INT32S sy, INT32S ex, INT32S ey);		//剪切图
/*-------------------------------- DRAW --------------------------------*/
void drawPoint		(INT32S x0, INT32S y0);								//画点
void drawLine		(INT32S x1, INT32S y1, INT32S x2, INT32S y2);		//画线
void drawCircle		(INT32S x0, INT32S y0, INT32S r);					//画圆
void drawEllipse	(INT32S x0, INT32S y0, INT32S rx, INT32S ry);		//画椭圆
void drawRectangle	(INT32S x1, INT32S y1, INT32S x2, INT32S y2);		//画矩形
void drawPolygon	(INT32S x[],INT32S y[],INT32S n);					//画多边形
void drawWave		(INT32S x[],INT32S y[],INT32S n);					//画曲线
void drawBezier		(INT32S x[],INT32S y[],INT32S n);					//画贝塞尔曲线
void drawGrid		(INT32S sx, INT32S sy, INT32S ex, INT32S ey, INT32S dx, INT32S dy);	//画网格
void drawCopy		(INT32S x0, INT32S y0, Mat<RGB>& gt);								//复制别的图
void fillRectangle	(INT32S sx, INT32S sy, INT32S ex, INT32S ey, ARGB color);			//填充单色
void fillFlood		(INT32S x0, INT32S y0, ARGB color);					//泛滥填充
void fillPolygon	(INT32S x[],INT32S y[],INT32S n);					//多边形填充
void drawChar		(INT32S x0, INT32S y0, char charac);				//显示字符
void drawString		(INT32S x0, INT32S y0, const char* str, INT32S n);	//显示字符串
void drawNum		(INT32S x0, INT32S y0, FP64 num);					//显示数字
/*-------------------------------- 二维变换 --------------------------------*/
void translate		(INT32S dx, INT32S dy);								//平移
void rotate			(FP64 theta,		INT32S x0, INT32S y0);			//旋转
void scale			(FP64 sx, FP64 sy,	INT32S x0, INT32S y0);			//缩放
```
![image](https://github.com/LiGuer/LiGu_Graphics/blob/master/example/LIGU-2D.jpg) 
![image](https://github.com/LiGuer/LiGu_Graphics/blob/master/example/LIGU-2D-Stm32.jpg) 

### <GraphicsND.h> 任意维图形学 + 科学制图:
```
/*---------------- 基础参数 ----------------*/
Graphics g;															//核心图形学类
Mat<Mat<int>> Z_Buffer;
Mat<> WindowSize{ 2,1 };											//窗口尺寸
static Mat<> TransformMat;											//变换矩阵
unsigned int FaceColor = 0xFFFFFF;
std::vector<Mat<>> LineSet, TriangleSet;
bool isLineTriangleSet = 0;
bool FACE = true, LINE = false;
/*---------------- 底层 ----------------*/
~GraphicsND() { ; }														//析构函数
GraphicsND(int width = 500, int height = 500, int Dim = 3) { init(width, height , Dim); }	//构造函数
void init (int width, int height, int Dim = 3);							//初始化
void clear(ARGB color);													//清屏
void value2pix	(int x0, int y0, int z0, int& x, int& y, int& z);		//点To像素 (<=3D)
void value2pix	(Mat<>& p0, Mat<int>& pAns);							//点To像素 (anyD)
bool setPix		(int x, int y, int z = 0, int size = -1);				//写像素 (正投影) (<=3D)
bool setPix		(Mat<int>& p0, int size = -1);							//写像素 (正投影) (anyD)
void setAxisLim	(Mat<>& pMin, Mat<>& pMax);								//设置坐标范围
void writeModel (const char* fileName);									//写模型文件
/*---------------- DRAW ----------------*/
// 0-D
void drawPoint		(double x0 = 0, double y0 = 0, double z0 = 0);	//画点 (<=3D)
void drawPoint		(Mat<>& p0);									//画点 (anyD)
// 1-D
void drawLine		(double sx0 = 0, double ex0 = 0, 
					 double sy0 = 0, double ey0 = 0, 
					 double sz0 = 0, double ez0 = 0);				//画直线 (<=3D)
void drawLine		(Mat<>& sp0, Mat<>& ep0);						//画直线 (anyD)
void drawPolyline	(Mat<>* p, int n, bool close = false);			//画折线
void drawPolyline	(Mat<>& y, double xmin, double xmax);			//画折线
void drawBezierLine	(Mat<> p[], int n);								//画Bezier曲线
// 2-D
void drawTriangle	(Mat<>& p1, Mat<>& p2, Mat<>& p3);						//画三角形
void drawTriangleSet(Mat<>& p1, Mat<>& p2, Mat<>& p3);						//画三角形集
void drawTriangleSet(Mat<>& p1, Mat<>& p2, Mat<>& p3, Mat<>&FaceVec);		//画三角形集
void drawRectangle	(Mat<>& sp, Mat<>& ep, Mat<>* direct = NULL);			//画矩形
void drawQuadrangle	(Mat<>& p1, Mat<>& p2, Mat<>& p3, Mat<>& p4);			//画四边形
void drawPolygon	(Mat<> p[], int n);										//画多边形
void drawCircle		(Mat<>& center, double r,				Mat<>* direct = NULL);		//画圆
void drawEllipse	(Mat<>& center, double rx, double ry,	Mat<>* direct = NULL);		//画椭圆
void drawSurface	(Mat<>& z, double xs, double xe, double ys, double ye, 
															Mat<>* direct = NULL);		//画曲面
void drawBezierFace	(Mat<> p[], int n);										//画贝塞尔曲面
// 3-D
void drawTetrahedron(Mat<>& p1, Mat<>& p2, Mat<>& p3, Mat<>& p4);		//画四面体
void drawCuboid		(Mat<>&pMin,Mat<>& pMax);							//画矩体
void drawPolyhedron	(Mat<>* p, int n);									//画多面体
void drawFrustum	(Mat<>& st, Mat<>& ed, double Rst, double Red, double delta = 36);	//画圆台
void drawCylinder	(Mat<>& st, Mat<>& ed, double r, double delta = 36);//画圆柱
void drawSphere		(Mat<>& center, double r, int delta = 36);			//画球
void drawSphere		(Mat<>& center, double r, double thetaSt, double thetaEd, 
						double phiSt = -PI / 2, double phiEd = PI / 2, int delta = 36);//画部分球
void drawSphere2	(Mat<>& center, double r, int n = 300);				//画球
void drawEllipsoid	(Mat<>& center, Mat<>& r);							//画椭球
void drawBody		(Mat<>& center, Mat<>& r);							//画曲体
void drawBezierBody	(Mat<> p[], int n);									//画Bezier曲体
void drawPipe		(Mat<>& st, Mat<>& ed, double Rst, double Red, int delta = 36);	//画平移体(粗细正多边形截面,线段路径)
void drawPipe		(Mat<>& st, Mat<>& ed, double R,               int delta = 36);	//画平移体(正多边形截面,线段路径)
void drawPipe		(Mat<>* p,  int N,	   double R,			   int delta = 36);	//画平移体(正多边形截面,any路径)
void drawPipe		(Mat<>& path,		   double R,			   int delta = 36);	//画平移体(正多边形截面,any路径)
void drawPipe		(Mat<>& st, Mat<>& ed, Mat<>& f);								//画平移体(any截面,线段路径)
void drawPipe		(Mat<>& path, Mat<>& f);										//画平移体(any截面,any路径)
void drawRotator	(Mat<>& zero, Mat<>& axis, Mat<>& f, int delta = 36, double st = 0, double ed = 2 * PI);	//画旋转体
void drawStairs		(Mat<>& zero, double Length, double Width, double Height, int Num);	//画阶梯
// Word
void drawChar	(Mat<>& p0, char charac);					//显示字符
void drawString	(Mat<>& p0, const char* str, int n);		//显示字符串
void drawNum	(Mat<>& p0, double num);					//显示数字
// any-D
void drawSuperLine	(Mat<>* p0);							//画线 any-D
void drawSuperCuboid(Mat<>& pMin, Mat<>& pMax);				//画立方体 any-D
void drawSuperSphere(Mat<>& center, double r);				//画球体 any-D
void drawGrid		(Mat<>& delta, Mat<>& max, Mat<>& min);	//画网格
// Other
void drawAxis(double Xmax = 0, double Ymax = 0, double Zmax = 0, bool negative = false);						//画坐标轴
void contour	(Mat<>& map, const int N);																		//画等高线
void contour	(Mat<>& map);																					//画等高面
void contour	(Mat<>& mapX, Mat<>& mapY, Mat<>& mapZ);
ARGB colorlist(double index, int model);																		//色谱
/*---------------- Transformation ----------------*/
static Mat<>& translate	(Mat<>& delta,										Mat<>& transMat = TransformMat);	//平移
static Mat<>& rotate	(double theta, Mat<>& center,						Mat<>& transMat = TransformMat);	//旋转 2D
static Mat<>& rotate	(Mat<>& rotateAxis, double theta, Mat<>& center,	Mat<>& transMat = TransformMat);	//旋转 3D
static Mat<>& rotate	(Mat<Mat<>>& rotateAxis, Mat<>& theta,Mat<>& center,Mat<>& transMat = TransformMat);	//旋转 4D
static Mat<>& scale		(Mat<>& ratio, Mat<>& center,						Mat<>& transMat = TransformMat);	//缩放
```  
![image](https://github.com/LiGuer/LiGu_Graphics/blob/master/example/高维空间/四维超立方图_四维超球.png) 
![image](https://github.com/LiGuer/LiGu_Graphics/blob/master/example/高维空间/四维超球-纬度.png) 
![image](https://github.com/LiGuer/LiGu_Graphics/blob/master/example/树.jpg) 
![image](https://github.com/LiGuer/LiGu_Graphics/blob/master/example/山海人.png)

### <RayTracing.h> 光线追踪:  

```
// 几何光学:
static Mat<>& reflect		(Mat<>& RayI, Mat<>& faceVec, Mat<>& RayO);								//反射
static Mat<>& refract		(Mat<>& RayI, Mat<>& faceVec, Mat<>& RayO, double rateI, double rateO);	//折射
static Mat<>& diffuseReflect(Mat<>& RayI, Mat<>& faceVec, Mat<>& RayO);								//漫反射
----------------------------------------------------------------------------------------
// 求交点:
double RayTriangle	(Mat<>& RaySt, Mat<>& Ray, Mat<>& p1, Mat<>& p2, Mat<>& p3);	//求交-射线与三角面
double RayPolygon	(Mat<>& RaySt, Mat<>& Ray, Mat<>* p,  int n);					//求交-射线与多边面
double RaySphere	(Mat<>& RaySt, Mat<>& Ray, Mat<>& center, double& R);			//求交-射线与球
double RayCuboid	(Mat<>& RaySt, Mat<>& Ray, Mat<>& p1, Mat<>& p2, Mat<>& p3);	//求交-射线与长方体
double RayCuboid	(Mat<>& RaySt, Mat<>& Ray, Mat<>& pmin, Mat<>& pmax);			//求交-射线与长方体 (轴对齐)
----------------------------------------------------------------------------------------
// 核心:
void paint(const char* fileName, int sampleSt = 0, int sampleEd = 0x7FFFFFFF);		//渲染
Mat<>& traceRay(Mat<>& RaySt, Mat<>& Ray, Mat<>& color, int level);			//追踪光线
double seekIntersection(Object& triangle, Mat<>& RaySt, Mat<>& Ray);		//求交点
```
![image](https://github.com/LiGuer/LiGu_Graphics/blob/master/example/RayTracing/RayTracing-Room-2.png) 
![image](https://github.com/LiGuer/LiGu_Graphics/blob/master/example/RayTracing/RayTracing-Room-3.png) 
![image](https://github.com/LiGuer/LiGu_Graphics/blob/master/example/RayTracing/色散火彩.png) 
![image](https://github.com/LiGuer/LiGu_Graphics/blob/master/example/RayTracing/RayTracing-Room-12885.png) 
![image](https://github.com/LiGuer/LiGu_Graphics/blob/master/example/RayTracing/RayTracing-Shpere.png) 


### <Fractal.h> 分形:
```
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
```

## Reference
* [1].Computer Graphics with OpenGL. Donald Hearn, M. Pauline Baker, Warren R. Carithers
* [2].正点原子-STM32