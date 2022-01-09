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

// Viewport stuff window 1
int windowW1 = 1400;
int windowH1 = 801;

double halfWindowW1 = windowW1 / 2.;
double halfWindowH1 = windowH1 / 2.;
double invW1 = 2. / windowW1;
double invH1 = 2. / windowH1;

// Viewport stuff window 2
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
    scale2 = 2 * scale / size_y;
    sinTheta = sin(theta);
    cosTheta = cos(theta);
}

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

    textureCoord.x = (pixelCoord.x * invW1 - 1) / viewScale1 + viewX1;
    textureCoord.y = ((windowH1 - pixelCoord.y) * invH1 - 1) / viewScale1 + viewY1;

    return textureCoord;
}

PixelCoord ttp(TextureCoord textureCoord) {
    PixelCoord pixelCoord;

    pixelCoord.x = ((textureCoord.x - viewX1) * viewScale1 + 1) * halfWindowH1;
    pixelCoord.y = ((textureCoord.y - viewY1) * viewScale1 + 1) * halfWindowW1;

    return pixelCoord;
}

// Util

bool operator==(const FractalCoord& p1, const FractalCoord& p2) {
    return (p1.x == p2.x && p1.y == p2.y);
}

bool operator!=(const FractalCoord& p1, const FractalCoord& p2) {
    return (p1.x != p2.x || p1.y != p2.y);
}