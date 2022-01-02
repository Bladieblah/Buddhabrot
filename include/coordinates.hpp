#ifndef COORDINATES_H
#define COORDINATES_H

// Fractal size
// #define size_x 5120
// #define size_y 2881
#define size_x 1400
#define size_y 801

extern double halfSize_x;
extern double halfSize_y;
extern double invHalfSize_x;
extern double invHalfSize_y;

// Fractal Positioning
extern double scale; // Half the y size
extern double dx;
extern double dy;
extern double theta;

extern double invScale;
extern double scale2;
extern double sinTheta;
extern double cosTheta;
extern double ratio_xy;
extern double ratio_yx;
extern double mutateSpread;

// Viewport stuff
extern int windowW;
extern int windowH;

extern double windowW2;
extern double windowH2;
extern double invW2;
extern double invH2;

// Texture positioning
extern double viewScale;
extern double viewX;
extern double viewY;
extern double drawScale;

extern double drawPower;


void updateFractalVars();

// Coordinate systems

// Coords in fractal array
typedef struct FractalCoord {
    int x, y;
} FractalCoord;

// Coords in mandelbrot space
typedef struct MandelCoord {
    double x, y;
} MandelCoord;

// Coords in texture space
typedef struct TextureCoord {
    double x, y;
} TextureCoord;

// Pixel coordinates on screen
typedef struct PixelCoord {
    int x, y;
} PixelCoord;

MandelCoord ftm(FractalCoord fractalCoord);
FractalCoord mtf(MandelCoord mandelCoord);
FractalCoord mtfMirror(MandelCoord mandelCoord);
TextureCoord mtt(MandelCoord mandelCoord);
MandelCoord ttm(TextureCoord textureCoord);
TextureCoord ptt(PixelCoord pixelCoord);
PixelCoord ttp(TextureCoord textureCoord);

bool operator==(const FractalCoord& p1, const FractalCoord& p2);
bool operator!=(const FractalCoord& p1, const FractalCoord& p2);

#endif