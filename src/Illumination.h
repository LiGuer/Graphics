#ifndef GRAPHICS_ILLMINATION_H
#define GRAPHICS_ILLMINATION_H

#include <vector>
#include <algorithm>

using namespace std;

/*
 *    Illumination Model
 */

namespace Illumination {

struct {
    vector<double>
        light = { -1 / sqrt(3), -1 / sqrt(3), 1 / sqrt(3) },
        view = { 0, 0, 1 };
    double 
        ka = 0.2, 
        kd = 0.9, 
        ks = 0.2, 
        alpha = 10;
} Phong_Parameter;

inline double Phong (double nx, double ny, double nz) {
    double a = sqrt(nx * nx + ny * ny + nz * nz);
    nx /= a;
    ny /= a;
    nz /= a;

    double 
        I = 0,
        b = nx * Phong_Parameter.light[0] +
            ny * Phong_Parameter.light[1] +
            nz * Phong_Parameter.light[2],
        rx = 2 * b * nx - Phong_Parameter.light[0],
        ry = 2 * b * ny - Phong_Parameter.light[1],
        rz = 2 * b * nz - Phong_Parameter.light[2];

    I = Phong_Parameter.ka +
        Phong_Parameter.kd * max(b, 0.0) +
        Phong_Parameter.ks * pow(max((
            rx * Phong_Parameter.view[0] +
            ry * Phong_Parameter.view[1] +
            rz * Phong_Parameter.view[2]
        ), 0.0), Phong_Parameter.alpha);

    return min(I, 1.0);
}

}

#endif