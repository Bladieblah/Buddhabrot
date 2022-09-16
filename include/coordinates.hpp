#ifndef COORDINATES_H
#define COORDINATES_H

// Fractal size
// #define size_x 5120
// #define size_y 2880

// Me
// #define size_x 6048
// #define size_y 3928

// Lun en ook ook Thek
// #define size_x 3840
// #define size_y 2160

// Thek
// #define size_x 3200
// #define size_y 1800
// Ook Thek
// #define size_x 3360
// #define size_y 2100

#define size_x 1400
#define size_y 801

#define size_x2 2600
#define size_y2 2000

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
extern double scaleDouble;
extern double sinTheta;
extern double cosTheta;
extern double ratio_xy;
extern double ratio_yx;
extern double mutateSpread;

// Fractal Positioning
extern double scale2; // Half the y size
extern double dx2;
extern double dy2;
extern double theta2;

extern double invScale2;
extern double scaleDouble2;
extern double sinTheta2;
extern double cosTheta2;
extern double ratio_xy2;
extern double ratio_yx2;

// Viewport stuff window 1
extern int windowW1;
extern int windowH1;

extern double halfWindowW1;
extern double halfWindowH1;
extern double invW1;
extern double invH1;

// Viewport stuff window 1
extern int windowW2;
extern int windowH2;

extern double halfWindowW2;
extern double halfWindowH2;
extern double invW2;
extern double invH2;

// Texture positioning window 1
extern double viewScale1;
extern double viewX1;
extern double viewY1;

// Texture positioning window 1
extern double viewScale2;
extern double viewX2;
extern double viewY2;

// Colormap stuff
extern double drawScale;
extern double drawPower;


void updateFractalVars();
void updateFractalVars2();

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
FractalCoord mtf2(MandelCoord mandelCoord);
FractalCoord mtfMirror(MandelCoord mandelCoord);
TextureCoord mtt(MandelCoord mandelCoord);
TextureCoord mtt2(MandelCoord mandelCoord);
MandelCoord ttm(TextureCoord textureCoord);
MandelCoord ttm2(TextureCoord textureCoord);
TextureCoord ptt(PixelCoord pixelCoord);
TextureCoord ptt2(PixelCoord pixelCoord);
PixelCoord ttp(TextureCoord textureCoord);

bool operator==(const FractalCoord& p1, const FractalCoord& p2);
bool operator!=(const FractalCoord& p1, const FractalCoord& p2);

MandelCoord operator+(const MandelCoord& m1, const MandelCoord& m2);
MandelCoord operator-(const MandelCoord& m1, const MandelCoord& m2);
MandelCoord operator/(const MandelCoord& m1, const MandelCoord& m2);
MandelCoord operator*(const MandelCoord& m1, const MandelCoord& m2);

MandelCoord operator+(const MandelCoord& m, const double& d);
MandelCoord operator-(const MandelCoord& m, const double& d);
MandelCoord operator*(const MandelCoord& m, const double& d);

MandelCoord mandelSin(MandelCoord m);
MandelCoord mandelCos(MandelCoord m);
MandelCoord mandelDiv(MandelCoord m1, MandelCoord m2);

#endif