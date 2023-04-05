#ifndef GRAPHICS_3D_H
#define GRAPHICS_3D_H

#include <queue>
#include <map>
#include <algorithm>
#include <functional>
#include "RGB.h"
#include "font.h"
#include "Graphics2D.h"
#include "Illumination.h"
#include "Transform.h"

namespace Graphics {

// Basic Geometry
bool drawPoint  (Mat<ARGB>& image, Mat<int>& Z_buf, int x = 0, int y = 0, int z = 0);
bool drawPoint  (Mat<ARGB>& image, Mat<int>& Z_buf, int x, int y, int z, double fx, double fy, double fz);
void drawLine   (Mat<ARGB>& image, Mat<int>& Z_buf, int sx = 0, int ex = 0, 
                                                    int sy = 0, int ey = 0, 
                                                    int sz = 0, int ez = 0);
void drawLine   (Mat<ARGB>& image, Mat<int>& Z_buf, vector<vector<int>>& p, bool close = false);
void drawTriangle 
                (Mat<ARGB>& image, Mat<int>& Z_buf, vector<int>& p1, vector<int>& p2, vector<int>& p3);
void drawTriangle
                (Mat<ARGB>& image, Mat<int>& Z_buf, vector<double>& p1, vector<double>& p2, vector<double>& p3);
void drawTriangleSet
                (Mat<ARGB>& image, Mat<int>& Z_buf, vector<vector<vector<int>>>& p);
void drawSphere (Mat<ARGB>& image, Mat<int>& Z_buf, int x, int y, int z, int r);
void drawFunction
                (Mat<ARGB>& image, Mat<int>& Z_buf, int xs, int ys, int zs, 
                                                    function<double(double, double, double)> f,
                                                    function<void(double, double, double, double&, double&, double&)> df);


}

#endif