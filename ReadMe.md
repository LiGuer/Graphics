# LiGu-Graphics 计算机图形学  
* <Graphics.h/cpp>          2D图形学类 (核心类)
* <Graphics3D.h/cpp>        3D图形学类 + 科学制图
* ------
* <ReadImg.exe>: 可实时动态显示所存储的图片。
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

### API <Graphics3D.h>:
```
/*---------------- 基础参数 ----------------*/
Graphics* g = NULL;														//核心图形学类
Mat<int> Z_index;
Mat<FP64> WindowSize{ 2,1 };											//窗口尺寸
static Mat<FP64> TransformMat;											//变换矩阵
/*---------------- 底层 ----------------*/
Graphics3D(int WindowSize_Width, int WindowSize_height);				//构造函数
~Graphics3D();															//析构函数
void init(int WindowSize_Width, int WindowSize_height);					//初始化
void value2pix(Mat<double>& p0, int& x, int& y, int& z);
/*---------------- DRAW ----------------*/
// 0-D
void drawPoint(double x0, double y0);									//画点
void drawPoint(Mat<double>& p0);										//画点
// 1-D
void drawLine(Mat<double>& sp, Mat<double>& ep);						//画直线
void drawPolyline(Mat<double>* p, int n);								//画折线
void contour(Mat<double>& map, const int N);							//画等高线
void contourface(Mat<double>& map, const int N);						//画等高线
// 2-D
void drawRectangle(Mat<double>& sp, Mat<double>& ep, Mat<double>* direct = NULL);//画矩形
void drawCircle(Mat<double>& center, double r, Mat<double>* direct = NULL);	//画圆
void drawEllipse(Mat<double>& center, double rx, double ry, Mat<double>* direct = NULL);//画椭圆
void drawPolygon(Mat<double> p[], int n);								//画多边形
void drawPlane(Mat<double>& p1, Mat<double>& p2, Mat<double>& p3);		//画平面
void drawSurface(Mat<double> z, double xs, double xe, double ys, double ye);//画曲面
void drawBezier();														//画贝塞尔曲面
// 3-D
void drawTetrahedron(Mat<double> p[]);									//画四面体
void drawCuboid(Mat<double> pMin, Mat<double> pMax);					//画矩体
void drawSphere(Mat<double>& center, double r);							//画球
void drawEllipsoid(Mat<double>& center, Mat<double>& r);				//画椭球
// Fill
void fill(Mat<double>& sp, Mat<double>& ep, RGB color);					//填充
void fillTriangle(Mat<double> p0[]);									//三角填充
void fillPolygon(Mat<double> p0[], int n);								//多边形填充
void fillflood(Mat<double>& p0, RGB color);								//泛滥填充
// Word
void drawChar(Mat<double>& p0, CHAR charac);							//显示字符
void drawString(Mat<double>& p0, const CHAR* str, INT32S n);			//显示字符串
void drawNum(Mat<double>& p0, FP64 num);								//显示数字
// Other
void grid();															//显示网格
RGB colorlist(double index, int model);									//色谱
/*---------------- Transformation ----------------*/
void translation(Mat<double>& delta, Mat<double>& transMat = TransformMat);	//平移
void rotate(Mat<double>& rotateAxis, double theta, Mat<double>& center, Mat<double>& transMat = TransformMat);//三维旋转-四元数
void scaling(Mat<double>& scale, Mat<double>& center, Mat<double>& transMat = TransformMat);//缩放
/*---------------- SET ----------------*/
void setAxisLim(Mat<double> pMin, Mat<double> pMax);					//判断坐标范围
```  
![image](https://github.com/LiGuer/LiGu_Graphics/blob/master/example/3D-surface.jpg) 
![image](https://github.com/LiGuer/LiGu_Graphics/blob/master/example/LIGU.jpg) 

## Reference
* [1].Computer Graphics with OpenGL. Donald Hearn, M. Pauline Baker, Warren R. Carithers
* [2].正点原子-STM32