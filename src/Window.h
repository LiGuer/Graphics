#ifndef WINDOW_H
#define WINDOW_H

#include "Graphics/GraphicsND.h"

/*
 * 视窗 : 三维图、三视图
 */
class WindowPix {
public:

    Mat<ARGB> imgXYZ;
    Mat<ARGB> imgXY;
    Mat<ARGB> imgXZ;
    Mat<ARGB> imgYZ;

    double scaleXYZ = 1, ratioXY;
    int size, sizeXY,
        OffsetX,
        OffsetY,
        OffsetZ;

    WindowPix() { ; }
    WindowPix(int _size, int _sizeXY) {
        init(_size, _sizeXY);
    }

    void init(int _size, int _sizeXY) {
        size = _size;
        sizeXY = _sizeXY;
        ratioXY = 1.0 * sizeXY / size;

        OffsetX = size / 2;
        OffsetY = size / 2;
        OffsetZ = size / 2;

        imgXYZ.zero(size, size);
        imgXY.zero(sizeXY, sizeXY);
        imgXZ.zero(sizeXY, sizeXY);
        imgYZ.zero(sizeXY, sizeXY);

        imgXYZ = 0xFF000000;
        imgXY  = 0xFF000000;
        imgXZ  = 0xFF000000;
        imgYZ  = 0xFF000000;
    }

    void operator=(ARGB color) {
        imgXYZ = color;
        imgXY = color;
        imgXZ = color;
        imgYZ = color;
    }

    inline void pixValue(double x, double y, double z, double& px, double& py, double& pz) {
        px = x * scaleXYZ + OffsetX;
        py = y * scaleXYZ + OffsetY;
        pz = z * scaleXYZ + OffsetZ;
    }

    void drawPoint(double x, double y, double z) {
        Graphics::drawPoint(imgXYZ,
            y * scaleXYZ + OffsetY,
            x * scaleXYZ + OffsetX,
            z * scaleXYZ + OffsetZ
        );
        Graphics::drawPoint(imgXY,
            (int)((y * scaleXYZ + OffsetY) * ratioXY),
            (int)((x * scaleXYZ + OffsetX) * ratioXY)
        );
        Graphics::drawPoint(imgXZ,
            (int)((z * scaleXYZ + OffsetZ) * ratioXY),
            (int)((x * scaleXYZ + OffsetX) * ratioXY)
        );
        Graphics::drawPoint(imgYZ,
            (int)((z * scaleXYZ + OffsetZ) * ratioXY),
            (int)((y * scaleXYZ + OffsetY) * ratioXY)
        );
    }

    inline void scale(double ratio) {
        scaleXYZ *= ratio;

        OffsetX *= (1 - ratio) * size / 2 / OffsetX + ratio;
        OffsetY *= (1 - ratio) * size / 2 / OffsetY + ratio;
        OffsetZ *= (1 - ratio) * size / 2 / OffsetZ + ratio;
    }
};


#endif