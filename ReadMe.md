# LiGu-Graphics 计算机图形学  
  
* 基础算法均已实现  
* 可用于教学演示  
* API:  
```
    /*---------------- 基础参数 ----------------*/  
	INT32U gSize[2] = { 2048,2048 };										//窗口尺寸  
	RGB* Map = NULL;														//图(底层)  
	RGB PaintColor;															//画笔颜色  
	INT32U PaintSize = 0;													//画笔大小  
	Mat<FP64> gM = { 3 };													//变换矩阵  
	/*---------------- 常数 ----------------*/  
	const RGB TRANSPARENT = 0xFFFFFFFF;										//RGB:透明  
	const INT32U FontSize = 16;												//字符大小  
	/*---------------- 底层 ----------------*/  
	~Graphics() { free(Map); }												//析构函数  
	void init();															//初始化  
	void clear();	 														//清屏  
	void clear(RGB color);	 												//清屏  
	void setPoint(INT32U x, INT32U y);										//底层画点  
	void setPoint(INT32U x, INT32U y, RGB color);							//底层画点(指定颜色)  
	RGB  readPoint(INT32U x, INT32U y); 									//读点   
	void PicWrite(const CHAR* filename);									//存图(底层)  
	/*---------------- DRAW ----------------*/  
	void drawPoint(INT32U x0, INT32U y0);									//画点  
	void drawLine(INT32U x1, INT32U y1, INT32U x2, INT32U y2);				//画线  
	void drawCircle(INT32U x0, INT32U y0, INT32U r);					    //画圆  
	void drawEllipse(INT32U x0, INT32U y0, INT32U rx,INT32U ry);			//画椭圆  
	void drawTriangle(INT32U x1, INT32U y1, INT32U x2, INT32U y2, INT32U x3, INT32U y3);//画三角形  
	void drawRectangle(INT32U x1, INT32U y1, INT32U x2, INT32U y2);		   	//画矩形  
	void drawWave(INT32U x[], INT32U y[], INT32U n);						//画曲线  
	void drawBezier(INT32S x[], INT32S y[], INT32U n);						//画贝塞尔曲线  
	void drawGrid(INT32U sx, INT32U sy, INT32U ex, INT32U ey, INT32U dx, INT32U dy);  
	void drawCopy(INT32U x0, INT32U y0, Graphics* gt);						//复制别的图  
	void fill(INT32U sx, INT32U sy, INT32U ex, INT32U ey);		   			//填充单色  
	void fill(INT32U sx, INT32U sy, INT32U ex, INT32U ey, RGB color);		//填充单色(指定颜色)  
	void floodfill(INT32U x0, INT32U y0, RGB color);						//泛滥填充  
	void drawChar(INT32U x0, INT32U y0, CHAR charac);						//显示字符  
	void drawString(INT32U x0, INT32U y0, const CHAR* str, INT32U n);		//显示字符串  
	void drawNum(INT32U x0, INT32U y0, FP64 num);							//显示数字  
	/*---------------- 二维变换 TRANSFORMATION ----------------*/  
	void translation(INT32S dx, INT32S dy);									//平移  
	void rotate(FP64 theta);												//旋转  
	void rotate(FP64 theta, INT32S x0, INT32S y0);							//旋转(基于基准点)  
	void scaling(FP64 sx, FP64 sy);											//缩放(>1直接完成变换)  
	void reflect();															//反射  
	void shear();															//错切  
	void confirmTrans();													//确认变换  
	/*---------------- SET ----------------*/  
	bool judgeOutRange(INT32U x0, INT32U y0);								//判断坐标是否过界  
	void setGSize(INT32U width, INT32U height);								//设置窗口尺寸  
```
![image](https://github.com/LiGuer/LiGu_Graphics/blob/master/LIGU.png) 
