#include <math.h>

#include "coordinates.hpp"

double halfSize_x = (double)size_x / 2;
double halfSize_y = (double)size_y / 2;
double invHalfSize_x = 1. / halfSize_x;
double invHalfSize_y = 1. / halfSize_y;

// Fractal Positioning
double scale = 1.3; // Half the y size
double dx = -0.5;
double dy = 0.;
double theta = 0.;

// double scale = 0.212578;
// double theta = -0.122461;
// double dx = -0.973908;
// double dy = 0.793633;

double invScale = 1. / scale;
double scale2 = 2 * scale / size_y;
double sinTheta = 0;
double cosTheta = 1;
double ratio_xy = (double)size_x / (double)size_y;
double ratio_yx = (double)size_y / (double)size_x;
double mutateSpread = 1e-3;

// Viewport stuff
int windowW = 1400;
int windowH = 801;

double windowW2 = windowW / 2.;
double windowH2 = windowH / 2.;
double invW2 = 2. / windowW;
double invH2 = 2. / windowH;

// Texture positioning
double viewScale = 1.;
double viewX = 0.;
double viewY = 0.;
double drawScale = 2.;

double drawPower = 1. / drawScale;

void updateFractalVars() {
    invScale = 1. / scale;
    scale2 = 2 * scale / size_y;
    sinTheta = sin(theta);
    cosTheta = cos(theta);
    mutateSpread = fmin(1e-3, 1e-2 * scale);
}

// Mandel <=> Fractal

// FractalCoord mtf(MandelCoord mandelCoord) {
//     FractalCoord fractalCoord;
//     double temp;

//     fractalCoord.x = mandelCoord.x - dx;
//     fractalCoord.y = mandelCoord.y - dy;

//     temp = (cosTheta * fractalCoord.x - sinTheta * fractalCoord.y) * invScale + 1;
//     fractalCoord.y = (cosTheta * fractalCoord.y + sinTheta * fractalCoord.x) * invScale;
//     fractalCoord.x = temp;

//     return fractalCoord;
// }

// MandelCoord ftm(FractalCoord fractalCoord) {
//     TextureCoord textureCoord;

//     textureCoord.x = fractalCoord.x * invHalfSize_x - 1;
//     textureCoord.y = fractalCoord.y * invHalfSize_y - 1;

//     return ttm(textureCoord);
// }

// MandelCoord ftm(FractalCoord fractalCoord) {
//     MandelCoord mandelCoord;

//     mandelCoord.x = fractalCoord.x * scale2 + dx - scale * ratio_xy;
//     mandelCoord.y = (size_y - fractalCoord.y) * scale2 + dy - scale;

//     return mandelCoord;
// }

FractalCoord mtf(MandelCoord mandelCoord) {
    FractalCoord fractalCoord;
    double tempx = mandelCoord.x - dx;
    double tempy = mandelCoord.y - dy;

    fractalCoord.x = (cosTheta * tempx - sinTheta * tempy + scale * ratio_xy) / scale2;
    fractalCoord.y = (cosTheta * tempy + sinTheta * tempx + scale) / scale2;

    return fractalCoord;
}

FractalCoord mtfMirror(MandelCoord mandelCoord) {
    FractalCoord fractalCoord;
    double tempx = mandelCoord.x - dx;
    double tempy = -mandelCoord.y - dy;

    fractalCoord.x = (cosTheta * tempx - sinTheta * tempy + scale * ratio_xy) / scale2;
    fractalCoord.y = (cosTheta * tempy + sinTheta * tempx + scale) / scale2;

    return fractalCoord;
}

// Mandel <=> Texture

TextureCoord mtt(MandelCoord mandelCoord) {
    TextureCoord textureCoord;
    double temp;

    textureCoord.x = mandelCoord.x - dx;
    textureCoord.y = mandelCoord.y - dy;

    temp = (cosTheta * textureCoord.x - sinTheta * textureCoord.y) * invScale;
    textureCoord.y = (cosTheta * textureCoord.y + sinTheta * textureCoord.x) * invScale;
    textureCoord.x = temp;

    return textureCoord;
}

MandelCoord ttm(TextureCoord textureCoord) {
    MandelCoord mandelCoord;
    double temp;

    mandelCoord.x = textureCoord.x * scale * ratio_xy;
    mandelCoord.y = textureCoord.y * scale;

    temp = (cosTheta * mandelCoord.x + sinTheta * mandelCoord.y) + dx;
    mandelCoord.y = (cosTheta * mandelCoord.y - sinTheta * mandelCoord.x) + dy;
    mandelCoord.x = temp;

    return mandelCoord;
}

// Texture <=> Pixel

TextureCoord ptt(PixelCoord pixelCoord) {
    TextureCoord textureCoord;

    textureCoord.x = (pixelCoord.x * invW2 - 1) / viewScale + viewX;
    textureCoord.y = ((windowH - pixelCoord.y) * invH2 - 1) / viewScale + viewY;

    return textureCoord;
}

PixelCoord ttp(TextureCoord textureCoord) {
    PixelCoord pixelCoord;

    pixelCoord.x = ((textureCoord.x - viewX) * viewScale + 1) * windowH2;
    pixelCoord.y = ((textureCoord.y - viewY) * viewScale + 1) * windowW2;

    return pixelCoord;
}

// Util

bool operator==(const FractalCoord& p1, const FractalCoord& p2) {
    return (p1.x == p2.x && p1.y == p2.y);
}

bool operator!=(const FractalCoord& p1, const FractalCoord& p2) {
    return (p1.x != p2.x || p1.y != p2.y);
}