# LiGu-Graphics 计算机图形学  
* <Graphics.h/cpp>          Pix图形学类 (核心类)
* <GraphicsND.h/cpp>        任意维图形学类 + 科学制图
* <DigitalImageProcessing.h>数字图像处理
* <ReadImg.exe>             实时动态显示图片
* ------[ example ]------ 
* <Boids.h>                 鸟群算法   

## API
### API <Graphics.h>:  
```
/*---------------- 基础参数 ----------------*/
INT32S gWidth = 100, gHeight = 100;										//窗口尺寸
RGBBASIC* Map = NULL;													//图
RGB PaintColor = 0xFFFFFF;												//画笔颜色
INT32S PaintSize = 0, FontSize=16;										//画笔大小//字符大小
Mat<FP64> TransMat{ 3 };												//变换矩阵
const RGB TRANSPARENT = 0xFFFFFFFF;										//RGB:透明
/*---------------- 底层 ----------------*/
Graphics(INT32S width, INT32S height) { init(width, height); }
~Graphics() { free(Map); }												//析构函数
void init();															//初始化
void init(INT32S width, INT32S height);									//初始化
void clear(RGB color);	 												//清屏
void setPoint(INT32S x, INT32S y, RGB color);							//底层画点
RGB  readPoint(INT32S x, INT32S y); 									//读点 
void PicWrite(const CHAR* filename);									//存图
bool judgeOutRange(INT32S x0, INT32S y0);								//判断过界
void transSelf();														//全图变换
void CutSelf(INT32S sx, INT32S sy, INT32S ex, INT32S ey);				//剪切图
/*---------------- DRAW ----------------*/
void drawPoint(INT32S x0, INT32S y0);									//画点
void drawLine(INT32S x1, INT32S y1, INT32S x2, INT32S y2);				//画线
void drawCircle(INT32S x0, INT32S y0, INT32S r);					    //画圆
void drawEllipse(INT32S x0, INT32S y0, INT32S rx,INT32S ry);			//画椭圆
void drawPolygon(INT32S x[], INT32S y[], INT32S n);						//画多边形
void drawRectangle(INT32S x1, INT32S y1, INT32S x2, INT32S y2);		   	//画矩形
void drawWave(INT32S x[], INT32S y[], INT32S n);						//画曲线
void drawBezier(INT32S x[], INT32S y[], INT32S n);						//画贝塞尔曲线
void drawGrid(INT32S sx, INT32S sy, INT32S ex, INT32S ey, INT32S dx, INT32S dy);//画网格
void drawCopy(INT32S x0, INT32S y0, RGBBASIC* gt, INT32S width, INT32S height);//复制别的图
void fillRectangle(INT32S sx, INT32S sy, INT32S ex, INT32S ey, RGB color);	//填充单色
void fillFlood(INT32S x0, INT32S y0, RGB color);						//泛滥填充
void fillPolygon(INT32S x[], INT32S y[], INT32S n);						//多边形填充
void drawChar(INT32S x0, INT32S y0, CHAR charac);						//显示字符
void drawString(INT32S x0, INT32S y0, const CHAR* str, INT32S n);		//显示字符串
void drawNum(INT32S x0, INT32S y0, FP64 num);							//显示数字
/*---------------- 二维变换 TRANSFORMATION ----------------*/
void translation(INT32S dx, INT32S dy);									//平移
void rotate(FP64 theta, INT32S x0, INT32S y0);							//旋转
void scaling(FP64 sx, FP64 sy, INT32S x0, INT32S y0);					//缩放
```
![image](https://github.com/LiGuer/LiGu_Graphics/blob/master/example/LIGU-2D.jpg) 
![image](https://github.com/LiGuer/LiGu_Graphics/blob/master/example/LIGU-2D-Stm32.jpg) 

### API <GraphicsND.h>:
```
/*---------------- 基础参数 ----------------*/
Graphics g;																//核心图形学类
Mat<int> Z_Buffer;
Mat<double> WindowSize{ 2,1 };											//窗口尺寸
static Mat<double> TransformMat;										//变换矩阵
unsigned int LineColor, FaceColor;
/*---------------- 底层 ----------------*/
GraphicsND() { ; }
GraphicsND(int width , int height) { init(width, height); }				//构造函数
~GraphicsND() { ; }														//析构函数
void init(int width, int height);										//初始化
void value2pix(int x0, int y0, int z0, int& x, int& y, int& z);			//点To像素 (<=3D)
void value2pix(Mat<double>& p0, Mat<int>& pAns);						//点To像素 (anyD)
bool setPix(int x, int y, int z = 0);									//写像素 (正投影)
void setAxisLim(Mat<double> pMin, Mat<double> pMax);					//设置坐标范围
/*---------------- DRAW ----------------*/
// 0-D
void drawPoint(double x0 = 0, double y0 = 0, double z0 = 0);			//画点 (<=3D)
void drawPoint(Mat<double>& p0);										//画点 (anyD)
// 1-D
void drawLine(double sx0 = 0, double ex0 = 0, double sy0 = 0, double ey0 = 0, double sz0 = 0, double ez0 = 0);					//画直线 (<=3D)
void drawLine(Mat<double>& sp0, Mat<double>& ep0);						//画直线 (anyD)
void drawPolyline(Mat<double>* p, int n, bool close = false);			//画折线
void drawBezierLine(Mat<double> p[], int n);							//画贝塞尔曲线
// 2-D
void drawTriangle(Mat<double>& p1, Mat<double>& p2, Mat<double>& p3, bool FACE = false, bool LINE = true);						//画三角形
void drawRectangle(Mat<double>& sp, Mat<double>& ep, Mat<double>* direct = NULL, bool FACE = false, bool LINE = true);			//画矩形
void drawQuadrilateral(Mat<double>& p1, Mat<double>& p2, Mat<double>& p3, Mat<double>& p4, bool FACE = false, bool LINE = true);//画四边形
void drawPolygon(Mat<double> p[], int n, bool FACE = false, bool LINE = true);													//画多边形
void drawCircle(Mat<double>& center, double r, Mat<double>* direct = NULL, bool FACE = false, bool LINE = true);				//画圆
void drawEllipse(Mat<double>& center, double rx, double ry, Mat<double>* direct = NULL);										//画椭圆
void drawSurface(Mat<double> z, double xs, double xe, double ys, double ye);													//画曲面
void drawBezierFace(Mat<double> p[], int n);																					//画贝塞尔曲面
// 3-D
void drawTetrahedron(Mat<double>& p1, Mat<double>& p2, Mat<double>& p3, Mat<double>& p4);										//画四面体
void drawCuboid(Mat<double>& pMin, Mat<double>& pMax);																			//画矩体
void drawPolyhedron(Mat<double>* p, int n, bool FACE = false, bool LINE = true);												//画多面体
void drawFrustum(Mat<double>& st, Mat<double>& ed, double Rst, double Red, double delta = 5, bool FACE = false, bool LINE = true);//画圆台
void drawCylinder(Mat<double>& st, Mat<double>& ed, double r, double delta = 5, bool FACE = false, bool LINE = true);			//画圆柱
void drawSphere(Mat<double>& center, double r, int delta = 5, bool FACE = false, bool LINE = true);								//画球
void drawSphere2(Mat<double>& center, double r, int n = 300);																	//画球
void drawEllipsoid(Mat<double>& center, Mat<double>& r);																		//画椭球
void drawBody(Mat<double>& center, Mat<double>& r);																				//画曲体
void drawBezierBody(Mat<double> p[], int n);																					//画贝塞尔曲体
// Word
void drawChar(Mat<double>& p0, char charac);							//显示字符
void drawString(Mat<double>& p0, const char* str, int n);				//显示字符串
void drawNum(Mat<double>& p0, double num);								//显示数字
// Other
void drawGrid();																												//画网格
void drawAxis(double Xmax = 0, double Ymax = 0, double Zmax = 0, bool negative = false);										//画坐标轴
void contour(Mat<double>& map, const int N);																					//画等高线
void contourface(Mat<double>& map, const int N);																				//画等高面
Graphics::ARGB colorlist(double index, int model);																				//色谱
/*---------------- Transformation ----------------*/
void translation(Mat<double>& delta, Mat<double>& transMat = TransformMat);														//平移
void rotate(Mat<double>& rotateAxis, double theta, Mat<double>& center, Mat<double>& transMat = TransformMat);					//三维旋转-四元数
void scaling(Mat<double>& scale, Mat<double>& center, Mat<double>& transMat = TransformMat);									//缩放
```  
![image](https://github.com/LiGuer/LiGu_Graphics/blob/master/example/LIGU-3D-surface.jpg) 
![image](https://github.com/LiGuer/LiGu_Graphics/blob/master/example/LIGU-3D-flower.jpg) 
![image](https://github.com/LiGuer/LiGu_Graphics/blob/master/example/LIGU-3D-Tree.png) 

## Reference
* [1].Computer Graphics with OpenGL. Donald Hearn, M. Pauline Baker, Warren R. Carithers
* [2].正点原子-STM32