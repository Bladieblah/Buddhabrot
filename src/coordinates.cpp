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
double scaleDouble = 2 * scale / size_y;
double sinTheta = 0;
double cosTheta = 1;
double ratio_xy = (double)size_x / (double)size_y;
double ratio_yx = (double)size_y / (double)size_x;

double scale2 = 1.3; // Half the y size
double dx2 = -0.5;
double dy2 = 0.;
double theta2 = 0.;

double invScale2 = 1. / scale2;
double scaleDouble2 = 2 * scale2 / size_y2;
double sinTheta2 = 0;
double cosTheta2 = 1;
double ratio_xy2 = (double)size_x2 / (double)size_y2;
double ratio_yx2 = (double)size_y2 / (double)size_x2;

// Viewport stuff window 1
// int windowW1 = 1400;
// int windowH1 = 801;
int windowW1 = 700;
int windowH1 = 401;

double halfWindowW1 = windowW1 / 2.;
double halfWindowH1 = windowH1 / 2.;
double invW1 = 2. / windowW1;
double invH1 = 2. / windowH1;

// Viewport stuff window 2
// int windowW2 = 1400;
// int windowH2 = 801;
int windowW2 = 700;
int windowH2 = 401;

double halfWindowW2 = windowW2 / 2.;
double halfWindowH2 = windowH2 / 2.;
double invW2 = 2. / windowW2;
double invH2 = 2. / windowH2;

// Texture positioning
double viewScale1 = 1.;
double viewX1 = 0.;
double viewY1 = 0.;

// Texture positioning
double viewScale2 = 1.;
double viewX2 = 0.;
double viewY2 = 0.;

// Rescale fractals before drawing with colormap
double drawScale = 2.;
double drawPower = 1. / drawScale;

void updateFractalVars() {
    invScale = 1. / scale;
    scaleDouble = 2 * scale / size_y;
    sinTheta = sin(theta);
    cosTheta = cos(theta);
}

void updateFractalVars2() {
    invScale2 = 1. / scale2;
    scaleDouble2 = 2 * scale2 / size_y2;
    sinTheta2 = sin(theta2);
    cosTheta2 = cos(theta2);
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

    fractalCoord.x = (cosTheta * tempx - sinTheta * tempy + scale * ratio_xy) / scaleDouble;
    fractalCoord.y = (cosTheta * tempy + sinTheta * tempx + scale) / scaleDouble;

    return fractalCoord;
}

FractalCoord mtf2(MandelCoord mandelCoord) {
    FractalCoord fractalCoord;
    double tempx = mandelCoord.x - dx2;
    double tempy = mandelCoord.y - dy2;

    fractalCoord.x = (cosTheta2 * tempx - sinTheta2 * tempy + scale2 * ratio_xy2) / scaleDouble2;
    fractalCoord.y = (cosTheta2 * tempy + sinTheta2 * tempx + scale2) / scaleDouble2;

    return fractalCoord;
}

MandelCoord ftm2(FractalCoord fractalCoord) {
    MandelCoord mandelCoord;
    double tempx = fractalCoord.x * scaleDouble2 - scale2 * ratio_xy2;
    double tempy = fractalCoord.y * scaleDouble2 - scale2;

    mandelCoord.x = (cosTheta2 * tempx + sinTheta2 * tempy) + dx2;
    mandelCoord.y = (cosTheta2 * tempy - sinTheta2 * tempx) + dy2;

    return mandelCoord;
}

FractalCoord mtfMirror(MandelCoord mandelCoord) {
    FractalCoord fractalCoord;
    double tempx = mandelCoord.x - dx;
    double tempy = -mandelCoord.y - dy;

    fractalCoord.x = (cosTheta * tempx - sinTheta * tempy + scale * ratio_xy) / scaleDouble;
    fractalCoord.y = (cosTheta * tempy + sinTheta * tempx + scale) / scaleDouble;

    return fractalCoord;
}

// Mandel <=> Texture

TextureCoord mtt(MandelCoord mandelCoord) {
    TextureCoord textureCoord;
    double tempx, tempy;

    tempx = mandelCoord.x - dx;
    tempy = mandelCoord.y - dy;

    textureCoord.x = (cosTheta * tempx - sinTheta * tempy) * invScale * ratio_yx;
    textureCoord.y = (cosTheta * tempy + sinTheta * tempx) * invScale;

    return textureCoord;
}

TextureCoord mtt2(MandelCoord mandelCoord) {
    TextureCoord textureCoord;
    double tempx, tempy;

    tempx = mandelCoord.x - dx2;
    tempy = mandelCoord.y - dy2;

    textureCoord.x = (cosTheta2 * tempx - sinTheta2 * tempy) * invScale2 * ratio_yx2;
    textureCoord.y = (cosTheta2 * tempy + sinTheta2 * tempx) * invScale2;

    return textureCoord;
}

MandelCoord ttm(TextureCoord textureCoord) {
    MandelCoord mandelCoord;
    double tempx, tempy;

    tempx = textureCoord.x * scale * ratio_xy;
    tempy = textureCoord.y * scale;

    mandelCoord.x = (cosTheta * tempx + sinTheta * tempy) + dx;
    mandelCoord.y = (cosTheta * tempy - sinTheta * tempx) + dy;

    return mandelCoord;
}

// Window 2
MandelCoord ttm2(TextureCoord textureCoord) {
    MandelCoord mandelCoord;
    double tempx, tempy;

    tempx = textureCoord.x * scale2 * ratio_xy2;
    tempy = textureCoord.y * scale2;

    mandelCoord.x = (cosTheta2 * tempx + sinTheta2 * tempy) + dx2;
    mandelCoord.y = (cosTheta2 * tempy - sinTheta2 * tempx) + dy2;

    return mandelCoord;
}

// Texture <=> Pixel

TextureCoord ptt(PixelCoord pixelCoord) {
    TextureCoord textureCoord;

    textureCoord.x = (pixelCoord.x * invW1 - 1) / viewScale1 + viewX1;
    textureCoord.y = ((windowH1 - pixelCoord.y) * invH1 - 1) / viewScale1 + viewY1;

    return textureCoord;
}

TextureCoord ptt2(PixelCoord pixelCoord) {
    TextureCoord textureCoord;

    textureCoord.x = (pixelCoord.x * invW2 - 1) / viewScale2 + viewX2;
    textureCoord.y = ((windowH2 - pixelCoord.y) * invH2 - 1) / viewScale2 + viewY2;

    return textureCoord;
}

PixelCoord ttp(TextureCoord textureCoord) {
    PixelCoord pixelCoord;

    pixelCoord.x = ((textureCoord.x - viewX1) * viewScale1 + 1) * halfWindowW1;
    pixelCoord.y = windowH1 - ((textureCoord.y - viewY1) * viewScale1 + 1) * halfWindowH1;

    return pixelCoord;
}

// Util

bool operator==(const FractalCoord& p1, const FractalCoord& p2) {
    return (p1.x == p2.x && p1.y == p2.y);
}

bool operator!=(const FractalCoord& p1, const FractalCoord& p2) {
    return (p1.x != p2.x || p1.y != p2.y);
}