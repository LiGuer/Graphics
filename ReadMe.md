# LiGu-Graphics 计算机图形学  
* 核心类: <Graphics.h/cpp>
* 科学制图类: <Plot.h/cpp>
* 3D图形学类: <Graphics3D.h/cpp>

## API:  
<Graphics.h>
```
/*---------------- 基础参数 ----------------*/
INT32S gWidth = 2048, gHeight=2048;										//窗口尺寸
RGB* Map = NULL;														//图
RGB PaintColor;															//画笔颜色
INT32S PaintSize = 0;													//画笔大小
Mat<FP64> gM;															//变换矩阵
/*---------------- 常数 ----------------*/
const RGB TRANSPARENT = 0xFFFFFFFF;										//RGB:透明
const INT32S FontSize = 16;												//字符大小
/*---------------- 底层 ----------------*/
~Graphics() { free(Map); }												//析构函数
void init();															//初始化
void clear(RGB color);	 												//清屏
void setPoint(INT32S x, INT32S y, RGB color);							//底层画点
RGB  readPoint(INT32S x, INT32S y); 									//读点 
void PicWrite(const CHAR* filename);									//存图
void confirmTrans();													//确认变换
/*---------------- DRAW ----------------*/
void drawPoint(INT32S x0, INT32S y0);									//画点
void drawLine(INT32S x1, INT32S y1, INT32S x2, INT32S y2);				//画线
void drawCircle(INT32S x0, INT32S y0, INT32S r);					    //画圆
void drawEllipse(INT32S x0, INT32S y0, INT32S rx,INT32S ry);			//画椭圆
void drawTriangle(INT32S x1, INT32S y1, INT32S x2, INT32S y2, INT32S x3, INT32S y3);//画三角形
void drawRectangle(INT32S x1, INT32S y1, INT32S x2, INT32S y2);		   	//画矩形
void drawWave(INT32S x[], INT32S y[], INT32S n);						//画曲线
void drawBezier(INT32S x[], INT32S y[], INT32S n);						//画贝塞尔曲线
void drawGrid(INT32S sx, INT32S sy, INT32S ex, INT32S ey, INT32S dx, INT32S dy);//画网格
void drawCopy(INT32S x0, INT32S y0, RGB* gt, INT32S width, INT32S height);//复制别的图
void fill(INT32S sx, INT32S sy, INT32S ex, INT32S ey, RGB color);		//填充单色
void floodfill(INT32S x0, INT32S y0, RGB color);						//泛滥填充
void drawChar(INT32S x0, INT32S y0, CHAR charac);						//显示字符
void drawString(INT32S x0, INT32S y0, const CHAR* str, INT32S n);		//显示字符串
void drawNum(INT32S x0, INT32S y0, FP64 num);							//显示数字
/*---------------- 二维变换 TRANSFORMATION ----------------*/
void translation(INT32S dx, INT32S dy);									//平移
void rotate(FP64 theta, INT32S x0, INT32S y0);							//旋转
void scaling(FP64 sx, FP64 sy);											//缩放(>1直接完成变换)
/*---------------- SET ----------------*/
bool judgeOutRange(INT32S x0, INT32S y0);								//判断坐标是否过界
void setSize(INT32S width, INT32S height);								//设置窗口尺寸
```
![image](https://github.com/LiGuer/LiGu_Graphics/blob/master/LIGU.png) 

## Reference
* [1].Computer Graphics with OpenGL. Donald Hearn, M. Pauline Baker, Warren R. Carithers
* [2].正点原子-STM32