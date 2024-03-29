#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <ctime>
#include <chrono>
#include <math.h>
#include <vector>
#include <thread>
#include <algorithm>
#include <iostream>

#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>

#include "colour.hpp"
#include "pcg.hpp"
#include "coordinates.hpp"

// For loading shader files
#define MAX_SOURCE_SIZE (0x100000)

#define iterations 20
#define particleCountSqrt 20
#define particleCount (particleCountSqrt * particleCountSqrt)
#define threadCount 6
uint64_t iterCount = 0;

int particleX, particleY;

size_t source_size;
char *source_str;

// Array to be drawn
#define uint32_max 4294967294
uint32_t data[size_y * size_x * 3];

// Colourmap stuff
uint32_t *colourMap;
int nColours = 10000;

// Kernel size for parallelisation
size_t global_item_size[2] = {(size_t)size_x, (size_t)size_y};
size_t local_item_size[2] = {(size_t)size_x, (size_t)size_y};

// Mouse
int mouse_x_down, mouse_x_up;
int mouse_y_down, mouse_y_up;
int mouse_state = GLUT_UP;

// Mandelbrot config
int kmax = 800;
double escape = 9;

bool transform = true;
bool withColormap = false;
bool showColorBar = true;
bool shouldDrawGrid = false;
bool targetedRandom = false;
bool naiveBuddha = true;
bool showDifference = false;
bool segmented = true;
bool showContrib = false;
bool updateTexture = true;
bool onlyBackground =  false;


#define thresholdCount 3
int32_t thresholds[thresholdCount] = {250, 5000, 32000};
double thresholdColours[thresholdCount * 3] = {
    // 0.4, 0., 0.4, 
    0.3, 0.0, 0.5,
    0.0, 0.5, 0.3,
    // 0.0, 0.25, 0.15,
    0.6, 0.5, 0.0,
    // 0.3, 0.25, 0.0,
    // 0.,0.,0.,
    // 1., 0., 0.
    // 0., 1., 0.,
    // 0., 0., 1.
};
uint64_t maxVals[thresholdCount] = {0};
uint64_t fractal[size_y * size_x * thresholdCount] = {0};
uint64_t frameFractal[size_y * size_x * thresholdCount] = {0};
uint64_t fractal2[size_y * size_x] = {0};
uint64_t frameFractal2[size_y * size_x] = {0};
uint64_t maxVal = 0;
uint64_t pixSum = 0;

double fractalContrib[size_y * size_x];
double maxContrib = 0;

std::chrono::high_resolution_clock::time_point frameTime = std::chrono::high_resolution_clock::now();

void echo(const char *text) {
    fprintf(stderr, "%s\n", text);
}

int viewIndex = 0;

double scaleHist[1000];
double thetaHist[1000];
double dxHist[1000];
double dyHist[1000];

MandelCoord **path;
FractalCoord **ppath;

double *xSamples;
double *ySamples;
int sampleSizeX = 6000;
int sampleSizeY = 5000;

// Points to track
typedef struct Particle {
    double x, y;
    uint64_t len, hitCount;
    double impact = 0;
} Particle;

Particle **particles;

bool pixelSort(const FractalCoord &p1, const FractalCoord &p2) {
    if (p1.x != p2.x) {
        return p1.x < p2.x;
    }
    return p1.y < p2.y;
}

void sortPixels(FractalCoord *pixels, int len) {
    std::sort(pixels, pixels + len, &pixelSort);
}

int countUnique(FractalCoord *pixels, int len) {
    int count = 1;
    sortPixels(pixels, len);

    for (int i=1; i<len; i++) {
        if (pixels[i-1] != pixels[i]) {
            count++;
        }
    }

    return count;
}

void mandelStep(MandelCoord *z, MandelCoord *c) {
    double temp = z->x * z->x - z->y * z->y + c->x;

    z->y = 2 * z->x * z->y + c->y;
    z->x = temp;
}

void mutatePoint(MandelCoord *z, double spread) {
    z->x += RANDN() * spread;
    z->y += RANDN() * spread;
}

Particle converge(double a, double b, int thread, int particleIndex) {
    MandelCoord z({0., 0.});
    MandelCoord c({a, b});
    MandelCoord oldZ({0., 0.});
    Particle result({0, 0, 0, 0});

    double c2 = a*a + b*b;
    int stepLimit = 2;
    
    // Check if c in main or secondary bulb
    if (256.0*c2*c2 - 96.0*c2 + 32.0*a - 3.0 < 0.0 || 16.0*(c2 + 2.0 * a + 1.0) - 1.0 < 0.0) {
        return result;
    }
    
    uint64_t j=0;

    int imax = thresholdCount;
    if (onlyBackground) {
        imax = 1;
    }

    for (int i=0; i<imax; i++) {
        for (; j<thresholds[i]; j++) {
            mandelStep(&z, &c);

            if (z.x * z.x + z.y * z.y > escape) {
                FractalCoord fc;
                int ind1, ind2;
                double weight = (uint64_t)pow(1 + log2(j + 1), 2);
                uint64_t hitCount = 0;
                uint64_t pixCount = 0;
                double impact = 0;

                for (int k=0; k<j; k++) {
                    fc = mtf(path[thread][k]);

                    if (fc.x < 0 || fc.x >= size_x || fc.y < 0 || fc.y >= size_y) {
                        continue;
                    }

                    ppath[thread][hitCount] = {fc.x, fc.y};
                    ind2 = size_x * fc.y + fc.x;

                    if (segmented) {
                        if ((int)(fc.x / particleX) * particleCountSqrt + (int)(fc.y / particleY) == particleIndex) {
                            hitCount++;
                            impact += weight / (1. + fractal2[ind2]);
                        }
                        else {
                            pixCount++;
                        }
                    }
                    else {
                        hitCount++;
                        impact += weight / (1. + fractal2[ind2]);
                    }

                    fractal2[ind2] += weight;

                    if (showDifference) {
                        frameFractal2[ind2] += weight;
                    }
                    
                    if (fractal2[ind2] > maxVal) {
                        maxVal = fractal2[ind2]; 
                    }

                    for (int l=i; l<=i; l++) {
                        ind1 = thresholdCount * ind2 + l;
                        // impact += 1. / (0.01 + fractal[ind1]);
                        fractal[ind1] += 1;
                        
                        if (showDifference) {
                            frameFractal[ind1] += 1;
                        }
                        
                        if (fractal[ind1] > maxVals[l]) {
                            maxVals[l] = fractal[ind1]; 
                        }
                    }
                }

                for (int k=0; k<j; k++) {
                    fc = mtfMirror(path[thread][k]);

                    if (fc.x < 0 || fc.x >= size_x || fc.y < 0 || fc.y >= size_y) {
                        continue;
                    }

                    ppath[thread][hitCount] = {fc.x, fc.y};
                    ind2 = size_x * fc.y + fc.x;
                    if (segmented) {
                        if ((int)(fc.x / particleX) * particleCountSqrt + (int)(fc.y / particleY) == particleIndex) {
                            hitCount++;
                            impact += weight / (1. + fractal2[ind2]);
                        }
                        else {
                            pixCount++;
                        }
                    }
                    else {
                        hitCount++;
                        impact += weight / (1. + fractal2[ind2]);
                    }

                    fractal2[ind2] += weight;


                    if (showDifference) {
                        frameFractal2[ind2] += weight;
                    }
                    
                    if (fractal2[ind2] > maxVal) {
                        maxVal = fractal2[ind2]; 
                    }

                    for (int l=i; l<=i; l++) {
                        ind1 = thresholdCount * ind2 + l;
                        // impact += 1. / (0.01 + fractal[ind1]);
                        fractal[ind1] += 1;

                        if (showDifference) {
                            frameFractal[ind1] += 1;
                        }
                        
                        if (fractal[ind1] > maxVals[l]) {
                            maxVals[l] = fractal[ind1]; 
                        }
                    }
                }

                pixCount += hitCount;
                pixSum += pixCount;

                int sourceInd = size_x * (int)(size_y * ((b + 1.3) / 2.6)) + (int)(size_x * ((a + 2.5) / 4.));
                if (sourceInd >= 0 && sourceInd < size_y * size_x) {
                    fractalContrib[sourceInd] += impact;

                    if (fractalContrib[sourceInd] > maxContrib) {
                        maxContrib = fractalContrib[sourceInd];
                    }
                }

                result = {a, b, j, hitCount, impact};

                return result;
            }

            if (z.x == oldZ.x && z.y == oldZ.y) {
                result.len = j;
                return result;
            }

            if (j == stepLimit) {
                oldZ = {z.x, z.y};
                stepLimit *= 2;
            }

            path[thread][j] = {z.x, z.y};
        }
    }


    return result;
}

void processFractal() {
    int i, j, k, l, ind1, ind2;

    double thr[thresholdCount];
    double floatval = 0;
    for (k=0; k<thresholdCount; k++) {
        thr[k] = 0.9 * (double)maxVals[k] + 1;
    }

    for (i=0; i<size_x; i++) {
        for (j=0; j<size_y; j++) {
            ind1 = thresholdCount * (size_x * j + i);
            ind2 = 3 * (size_x * j + i);

            if (!showDifference) {
                for (l=0; l<3; l++) {
                    floatval = sqrt(fractal[ind1] / thr[0]) * thresholdColours[l];
                    // data[ind2 + l] = (uint64_t)(uint32_max * fmin(1., sqrt(fractal[ind1] / thr[0])) * thresholdColours[l]);
                    for (k=1; k<thresholdCount; k++) {
                        floatval += sqrt(fractal[ind1 + k] / thr[k]) * thresholdColours[3 * k + l];
                        // data[ind2 + l] += (uint64_t)(uint32_max * fmin(1., sqrt(fractal[ind1 + k] / thr[k])) * thresholdColours[3 * k + l]);
                    }

                    floatval = fmin(1., floatval);
                    data[ind2 + l] = (uint64_t)(uint32_max * floatval);
                }
            }
            else {
                for (l=0; l<3; l++) {
                    data[ind2 + l] = (uint64_t)(uint32_max * fmin(1., sqrt(frameFractal[ind1] / thr[0])) * thresholdColours[l]);
                    for (k=1; k<thresholdCount; k++) {
                        data[ind2 + l] += (uint64_t)(uint32_max * fmin(1., sqrt(frameFractal[ind1 + k] / thr[k])) * thresholdColours[3 * k + l]);
                    }
                }
            }
        }
    }
}

void processFractal2() {
    int i, j, ind, colInd;

    double thr = (double)maxVal;

    for (i=0; i<size_x; i++) {
        for (j=0; j<size_y; j++) {
            ind = size_x * j + i;

            if (showColorBar && i < 70) {
                ind *= 3;
                colInd = 3 * (int)(j / (double)size_y * nColours);
                data[ind + 0] = colourMap[colInd + 0];
                data[ind + 1] = colourMap[colInd + 1];
                data[ind + 2] = colourMap[colInd + 2];
                continue;
            }

            if (!showDifference) {
                colInd = 3 * (int)(nColours * pow((double)fractal2[ind] / thr, drawPower));
            }
            else {
                colInd = 3 * (int)(nColours * pow((double)frameFractal2[ind] / thr, drawPower));
            }

            ind *= 3;
            for (int k=0; k<3; k++) {
                data[ind + k] = colourMap[colInd + k];
            }
        }
    }
}

void processContrib() {
    int i, j, ind, colInd;

    double thr = (double)maxContrib;

    for (i=0; i<size_x; i++) {
        for (j=0; j<size_y; j++) {
            ind = size_x * j + i;

            if (showColorBar && i < 70) {
                ind *= 3;
                colInd = 3 * (int)(j / (double)size_y * nColours);
                data[ind + 0] = colourMap[colInd + 0];
                data[ind + 1] = colourMap[colInd + 1];
                data[ind + 2] = colourMap[colInd + 2];
                continue;
            }

            if (!showDifference) {
                colInd = 3 * (int)(nColours * pow((double)fractalContrib[ind] / thr, drawPower));
            }
            else {
                colInd = 3 * (int)(nColours * pow((double)fractalContrib[ind] / thr, drawPower));
            }

            ind *= 3;
            for (int k=0; k<3; k++) {
                data[ind + k] = colourMap[colInd + k];
            }
        }
    }
}

MandelCoord mutateParticle(int thread, int particle) {
    MandelCoord result;

    if (UNI() < 0.5) {
        result.x = particles[thread][particle].x + RANDN() / (1 + pow(particles[thread][particle].len, 1.2));
        result.y = particles[thread][particle].y + RANDN() / (1 + pow(particles[thread][particle].len, 1.2));
        // result.x = 2 * particles[thread][particle].x + RANDN() * fmin(scale, 1e-1);
        // result.y = 2 * particles[thread][particle].y + RANDN() * fmin(scale, 1e-1);
    }
    else {
        // result.x = 4.5 * UNI() - 2.6;
        // result.y = 3.0 * UNI() - 1.5;
        int si = UNI() * sampleSizeX;
        int sj = UNI() * sampleSizeY;
        double r, rx, ry;

        if (si == 0)
            rx = fabs(xSamples[sampleSizeX * sj + (si + 1)] - xSamples[sampleSizeX * sj + si]);
        else
            rx = fabs(xSamples[sampleSizeX * sj + (si - 1)] - xSamples[sampleSizeX * sj + si]);

        if (sj == 0)
            ry = fabs(xSamples[sampleSizeX * (sj + 1) + si] - xSamples[sampleSizeX * sj + si]);
        else
            ry = fabs(xSamples[sampleSizeX * (sj - 1) + si] - xSamples[sampleSizeX * sj + si]);

        r = fmax(rx, ry);
        if (r == 0) {
            r = fmax(rx, ry);
        }
        result.x = xSamples[sampleSizeX * sj + si] + RANDN() * r;
        result.y = ySamples[sampleSizeX * sj + si] + RANDN() * r;
    }

    return result;
}

void _metrobrot(int thread) {
    int i, j;
    double prob;
    Particle result;

    for (i=0; i<particleCount; i++) {
        for (j=0; j<iterations; j++) {
            MandelCoord testPoint = mutateParticle(thread, i);
            result = converge(testPoint.x, testPoint.y, thread, i);

            if (result.hitCount == 0) {
                continue;
            }

            particles[thread][i].impact *= 0.95;
            prob = pow(result.hitCount / (double)particles[thread][i].hitCount, 4)
                * pow(result.len / (double)particles[thread][i].len, 1)
                * pow(result.impact / (double)particles[thread][i].impact, 2)
            ;

            if (UNI() < prob) {
                particles[thread][i] = result;
            }
        }
    }
}

void _mandelbrot(int thread) {
    int si, sj;
    double a, b;
    double prob = -1;
    Particle result;
    bool _segmented = segmented;
    segmented = false;

	for (int i=0; i<iterations * particleCount; i++) {
        si = UNI() * sampleSizeX;
        sj = UNI() * sampleSizeY;

        a = xSamples[sampleSizeX * sj + si];
        b = ySamples[sampleSizeX * sj + si];
        // if (UNI() > prob) {
        //     a = 4.5 * UNI() - 2.6;
        //     b = 3.0 * UNI() - 1.5;
        // }
        // else {
        //     a += RANDN() * mutateSpread;
        //     b += RANDN() * mutateSpread;
        // }
        result = converge(a, b, thread, 0);
        // if (result.hitCount > 0) {
        //     prob = sqrt(result.len / (double)thresholds[thresholdCount-1] * sqrt(1 - countUnique(ppath[thread], result.hitCount) / (double)result.hitCount) + result.impact);
        // }
        // else {
        //     prob = 0;
        // }
    }

    segmented = _segmented;
}

void mandelbrot() {
	std::thread *tt = new std::thread[threadCount-1];
	
    if (naiveBuddha) {
        for (int i=0; i<threadCount-1; i++) {
            tt[i] = std::thread(_mandelbrot, i);
        }
        
        _mandelbrot(threadCount - 1);
    }
    else {
        for (int i=0; i<threadCount-1; i++) {
            tt[i] = std::thread(_metrobrot, i);
        }
        
        _metrobrot(threadCount - 1);
    }
	
	for (int i = 0; i < threadCount-1; ++i) {
		tt[i].join();
	}
	
	delete [] tt;
}

void _initParticles(int thread) {
    int j, k, si, sj;
    double a, b;
    Particle result;

    for (j=0; j<particleCount; j++) {
        for (k=0; k<10000; k++) {
            // a = 4.5 * UNI() - 2.6;
            // b = 3.0 * UNI() - 1.5;
            si = UNI() * sampleSizeX;
            sj = UNI() * sampleSizeY;

            a = xSamples[sampleSizeX * sj + si];
            b = ySamples[sampleSizeX * sj + si];

            result = converge(a, b, thread, j);

            if (result.hitCount > 0) {
                particles[thread][j] = result;
                break;
            }
        }

        if (k == 100000) {
            echo("Well ffs >:U");
            particles[thread][j] = result;
        }
    }
}

void initParticles() {
	std::thread *tt = new std::thread[threadCount-1];
	
    for (int i=0; i<threadCount-1; i++) {
        tt[i] = std::thread(_initParticles, i);
    }
    
    _initParticles(threadCount - 1);
	
	for (int i = 0; i < threadCount-1; ++i) {
		tt[i].join();
	}
	
	delete [] tt;
}

double computeDistance(int pathLength) {
    double dist = 1e10;
    double test;

    for (int i=0; i<pathLength; i++) {
        test = sqrt(pow(dx - path[0][i].x, 2) + pow(dy - path[0][i].y, 2));

        if (test < dist) {
            dist = test;
        }
    }

    return dist;
}

// Functions

void makeColourmap() {
    std::vector<float> x = {0., 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.6, 0.7, 1.};
    std::vector< std::vector<float> > y = {
        {0,0,0}, // #000000 Black
        {26,17,36}, // #1a1124 Darker purple
        {52,21,56}, // #341538 Dark purple
        {33,130,133}, // #218285 Teal
        {47,189,133}, // #2fbd85 Green
        {26,17,36}, // #1a1124 Darker purple
        {200,40,187}, // #c828bb Bright pink
        {241,249,244}, // #f1f9f4 White
        {210,245,255}, // #d2f5ff Dyme Blue
        {26,17,36} // #1a1124 Darker Purple
    };
    
    // std::vector<float> x = {0., 0.3333333, 0.6666666, 1.};
    // std::vector< std::vector<float> > y = {
    //     {0, 0, 64},
    //     {0, 255, 192},
    //     {64, 192, 0},
    //     {0, 0, 64}
    // };
    
    // std::vector<float> x = {0., 0.5, 1.};
    // std::vector< std::vector<float> > y = {
    //     {0, 0, 0},
    //     {0, 255, 255},
    //     {255, 255, 255}
    // };

    Colour col(x, y, nColours);
    colourMap = (uint32_t *)malloc(3 * nColours * sizeof(uint32_t));
    col.apply(colourMap);
    
    // Write colourmap to GPU
    // ret = clEnqueueWriteBuffer(command_queue, mapmobj, CL_TRUE, 0, 3*nColours*sizeof(float), colourMap, 0, NULL, NULL);
}

void prepare() {
    pcg32_srandom(time(NULL) ^ (intptr_t)&printf, (intptr_t)&kmax); // Seed pcg
    
    path = (MandelCoord **)malloc(threadCount * sizeof(MandelCoord *));
    ppath = (FractalCoord **)malloc(threadCount * sizeof(FractalCoord *));

    particles = (Particle **)malloc(threadCount * sizeof(Particle *));

    for (int i=0; i<threadCount; i++) {
        path[i] = (MandelCoord *)malloc(thresholds[thresholdCount-1] * sizeof(MandelCoord));
        ppath[i] = (FractalCoord *)malloc(2 * thresholds[thresholdCount-1] * sizeof(FractalCoord));
        particles[i] = (Particle *)malloc(particleCount * sizeof(Particle));
    }

    srand(time(NULL));
}

void cleanup() {
    free(colourMap);
    // free(fractal);
    free(source_str);
    free(path);
    free(ppath);

    free(xSamples);
    free(ySamples);
}

void step() {
    if (showDifference) {
        for (int i=0; i<size_x*size_y; i++) {
            frameFractal2[i] = 0;
            for (int j=0; j<thresholdCount; j++) {
                frameFractal[thresholdCount * i + j] = 0;
            }
        }
    }

    mandelbrot();

    if (showContrib) {
        processContrib();
    }
    else 
    if (withColormap) {
        processFractal2();
    }
    else {
        processFractal();
    }
}

void drawBox() {
    double xc = 2 * mouse_x_down / (double)windowW - 1;
    double yc = 1 - 2 * mouse_y_down / (double)windowH;

    double dx1 = 2 * (mouse_x_up - mouse_x_down) / (double)windowW;
    double dy1 = -2 * (mouse_y_up - mouse_y_down) / (double)windowH;

    double dx2 = 2 * (mouse_y_up - mouse_y_down) * ratio_xy / windowW;
    double dy2 = 2 * (mouse_x_up - mouse_x_down) * ratio_xy / windowH;

    glColor4f(1,1,1,1);

    glBegin(GL_LINE_STRIP);
        glVertex2f(xc + dx1 + dx2, yc + dy1 + dy2);
        glVertex2f(xc + dx1 - dx2, yc + dy1 - dy2);
        glVertex2f(xc - dx1 - dx2, yc - dy1 - dy2);
        glVertex2f(xc - dx1 + dx2, yc - dy1 + dy2);
        glVertex2f(xc + dx1 + dx2, yc + dy1 + dy2);
    glEnd();
}

void drawGrid() {
    glBegin(GL_LINES);
        glVertex2f(-1,0); glVertex2f(1,0);
        glVertex2f(-1,0.5); glVertex2f(1,0.5);
        glVertex2f(-1,-0.5); glVertex2f(1,-0.5);
        glVertex2f(0,-1); glVertex2f(0,1);
        glVertex2f(0.5,-1); glVertex2f(0.5,1);
        glVertex2f(-0.5,-1); glVertex2f(-0.5,1);
    glEnd();
}

void display() {
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    step();
    iterCount += iterations * threadCount;
    
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span1 = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - frameTime);
    std::chrono::duration<double> time_span2 = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    
    double avgImpact = 0;
    for (int i=0; i<threadCount; i++) {
        for (int j=0; j<particleCount; j++) {
            avgImpact += particles[i][j].impact;
        }
    }
    avgImpact /= threadCount * particleCount;

    fprintf(stderr, "\rIterations = %llu, pixSum = %llu, Frame time = %.4g, Step time = %.4g, impact = %g, maxContrib = %f", 
        iterCount, pixSum, time_span1.count(), time_span2.count(), avgImpact, maxContrib);
    
    frameTime = std::chrono::high_resolution_clock::now();

    if (updateTexture) {
        glClearColor( 0, 0, 0, 1 );
        glColor3f(1, 1, 1);
        glClear( GL_COLOR_BUFFER_BIT );

        glEnable (GL_TEXTURE_2D);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

        glTexImage2D (
            GL_TEXTURE_2D,
            0,
            GL_RGB,
            size_x,
            size_y,
            0,
            GL_RGB,
            GL_UNSIGNED_INT,
            &data[0]
        );

        glPushMatrix();
        glScalef(viewScale, viewScale, 1.);
        glTranslatef(-viewX, -viewY, 0.);

        glBegin(GL_QUADS);
            glTexCoord2f(0.0f, 0.0f); glVertex2f(-1.0, -1.0);
            glTexCoord2f(1.0f, 0.0f); glVertex2f( 1.0, -1.0);
            glTexCoord2f(1.0f, 1.0f); glVertex2f( 1.0,  1.0);
            glTexCoord2f(0.0f, 1.0f); glVertex2f(-1.0,  1.0);
        glEnd();

        glDisable (GL_TEXTURE_2D);
        glPopMatrix();

        if (shouldDrawGrid) {
            drawGrid();
        }

        if (mouse_state == GLUT_DOWN) {
            drawBox();
        }
        glFlush();
        glutSwapBuffers();
    }
}

void writeData() {
    int i, j, k;
    FILE *outFile;
    char filename[100];

    echo("Writing filename");
    sprintf(filename, "data/weighted_%d_%d_%.6f_%.6f_%.6f_%.6f.csv", size_x, size_y, scale, theta, dx, dy);
    
    echo("Opening File");
    outFile = fopen(filename, "w");

    if (outFile == NULL) {
        fprintf(stderr, "\nError opening file\n");
        return;
    }

    echo("Writing to file");
    k = 0;
    for (j=0; j<size_y; j++) {
        fprintf(outFile, "%llu", fractal2[k]);
        k++;

        for (i=1; i<size_x; i++) {
            fprintf(outFile, ",%llu", fractal2[k]);
            k++;
        }
        fprintf(outFile, "\n");
    }

    fclose(outFile);

    for (int l=0; l<thresholdCount; l++) {
        fprintf(stderr, "Writing filename %d\n", l);
        sprintf(filename, "data/raw_%d_%d_%.6f_%.6f_%.6f_%.6f_%d.csv", size_x, size_y, scale, theta, dx, dy, thresholds[l]);
        
        echo("Opening File");
        outFile = fopen(filename, "w");

        if (outFile == NULL) {
            fprintf(stderr, "\nError opening file\n");
            return;
        }

        echo("Writing to file");
        k = 0;
        for (j=0; j<size_y; j++) {
            fprintf(outFile, "%llu", fractal[thresholdCount * k + l]);
            k++;

            for (i=1; i<size_x; i++) {
                fprintf(outFile, ",%llu", fractal[thresholdCount * k + l]);
                k++;
            }
            fprintf(outFile, "\n");
        }

        fclose(outFile);
    }
}

void readSamples() {
    int i, j, k;
    FILE *inFile;
    char filename[100];

    xSamples = (double *)malloc(sampleSizeX * sampleSizeY * sizeof(double));
    ySamples = (double *)malloc(sampleSizeX * sampleSizeY * sizeof(double));

    sprintf(filename, "XX.csv");
    inFile = fopen(filename, "r");

    if (inFile == NULL) {
        fprintf(stderr, "\nError opening XX file\n");
        return;
    }

    k = 0;
    for (j=0; j<sampleSizeY; j++) {
        for (i=0; i<sampleSizeX-1; i++) {
            fscanf(inFile, "%lf,", &xSamples[k]);
            k++;
        }
        fscanf(inFile, "%lf\n", &xSamples[k]);
        k++;
    }

    fclose(inFile);

    sprintf(filename, "YY.csv");
    inFile = fopen(filename, "r");

    if (inFile == NULL) {
        fprintf(stderr, "\nError opening YY file\n");
        return;
    }

    k = 0;
    for (j=0; j<sampleSizeY; j++) {
        for (i=0; i<sampleSizeX-1; i++) {
            fscanf(inFile, "%lf,", &ySamples[k]);
            k++;
        }
        fscanf(inFile, "%lf\n", &ySamples[k]);
        k++;
    }

    fclose(inFile);
}

void readData() {
    int i, j, k;
    FILE *inFile;
    char filename[100];

    sprintf(filename, "data/weighted_%d_%d.csv", size_x, size_y);
    inFile = fopen(filename, "r");

    if (inFile == NULL) {
        fprintf(stderr, "\nError opening file\n");
        return;
    }

    k = 0;
    for (j=0; j<size_y; j++) {
        for (i=0; i<size_x-1; i++) {
            fscanf(inFile, "%llu,", &fractal2[k]);
            k++;
        }
        fscanf(inFile, "%llu\n", &fractal2[k]);
        k++;
    }

    fclose(inFile);
}

void clearData() {
    int i,j;

    for (i=0; i<size_x * size_y; i++) {
        fractal2[i] = 0;
        fractalContrib[i] = 0;
        for (j=0; j<thresholdCount; j++) {
            fractal[thresholdCount * i + j] = 0;
        }
    }


    maxVal = 0;
    pixSum = 0;
    iterCount = 0;
    for (j=0; j<thresholdCount; j++) {
        maxVals[j] = 0;
    }
}

void setCoordinates(double _scale, double _theta, double _dx, double _dy) {
    scaleHist[viewIndex] = scale;
    thetaHist[viewIndex] = theta;
    dxHist[viewIndex] = dx;
    dyHist[viewIndex] = dy;

    fprintf(stderr, "\nscale (%f -> %f)\n", scale, _scale);
    fprintf(stderr, "theta (%f -> %f)\n", theta, _theta);
    fprintf(stderr, "dx (%f -> %f)\n", dx, _dx);
    fprintf(stderr, "dy (%f -> %f)\n", dy, _dy);

    double dx1 = _scale / size_y * size_x;
    double dx2 = scale / size_y * size_x;
    fprintf(stderr, "x (%f - %f) -> (%f - %f)\n", dx - dx2, dx + dx2, _dx - dx1, _dx + dx1);

    scale = _scale;
    theta = _theta;
    dx = _dx;
    dy = _dy;

    viewIndex++;
    
    echo("Initialising particles");
    initParticles();
    echo("Done");

    updateFractalVars();

    clearData();
    initParticles();
}

void selectRegion() {
    if (!transform || mouse_state != GLUT_DOWN) {
        return;
    }

    PixelCoord pCenter({mouse_x_down, mouse_y_down});
    PixelCoord pTop({mouse_x_up, mouse_y_up});

    MandelCoord mCenter = ttm(ptt(pCenter));
    MandelCoord mTop = ttm(ptt(pTop));

    setCoordinates(
        sqrt(pow(mTop.x - mCenter.x, 2) + pow(mTop.y - mCenter.y, 2)),
        atan2(mTop.x - mCenter.x, mTop.y - mCenter.y),
        mCenter.x,
        mCenter.y
    );
}

void revertCoordinates() {
    if (viewIndex == 0) {
        return;
    }

    viewIndex--;
    
    scale = scaleHist[viewIndex];
    theta = thetaHist[viewIndex];
    dx = dxHist[viewIndex];
    dy = dyHist[viewIndex];

    updateFractalVars();
    clearData();
    
    if (viewIndex > 0) {
        echo("\nInitialising particles");
        initParticles();
        echo("Done");
    }
}


void keyPressed(unsigned char key, int x, int y) {
    switch (key) {
        case 'w':
            viewScale *= 1.5;
            step();
            break;
        case 's':
            viewScale /= 1.5;
            step();
            break;
        case 'p':
            glutIdleFunc(&display);
            break;
        case 'e':
            // echo("Pressing e");
            glutPostRedisplay();
            break;
        case 't':
            transform = !transform;
            fprintf(stderr, "Set transform to %d\n", transform);
            break;
        case 'c':
            withColormap = !withColormap;
            break;
        case 'b':
            showColorBar = !showColorBar;
            break;
        case 'B':
            onlyBackground = !onlyBackground;
            break;
        case 'g':
            shouldDrawGrid = !shouldDrawGrid;
            break;
        case 'o':
            fprintf(stderr, "\nmaxVal = %llu, pixSum = %llu", maxVal, pixSum);
            for (int k=0; k<thresholdCount; k++){
                fprintf(stderr, ", maxval %d = %llu", k, maxVals[k]);
            }
            fprintf(stderr, "\n");
            break;
        case 'm':
            drawScale *= 1.1;
            drawPower = 1. / drawScale;
            break;
        case 'n':
            drawScale /= 1.1;
            drawPower = 1. / drawScale;
            break;
        case 'D':
            writeData();
            break;
        case 'a':
            selectRegion();
            break;
        case 'r':
            viewX = 0.;
            viewY = 0.;
            viewScale = 1.;
            drawScale = 2.;
            drawPower = 1. / drawScale;
            break;
        case 'i':
            initParticles();
            break;
        case 'z':
            revertCoordinates();
            break;
        case 'j':
            targetedRandom = !targetedRandom;
            fprintf(stderr, "Set targetedRandom to %d\n", targetedRandom);
            // clearData();
            break;
        case 'x':
            segmented = !segmented;
            fprintf(stderr, "Set segmented to %d\n", segmented);
            break;
        case 'l':
            naiveBuddha = !naiveBuddha;
            fprintf(stderr, "Set naiveBuddha to %d\n", naiveBuddha);
            // clearData();
            break;
        case 'd':
            showDifference = !showDifference;
            break;
        case 'u':
            showContrib = !showContrib;
            break;
        case 'U':
            updateTexture = !updateTexture;
            break;
        case 'q':
        	cleanup();
        	fprintf(stderr, "\n");
            exit(0);
            break;
        default:
            break;
    }
}

void specialKeyPressed(int key, int x, int y) {
    switch (key) {
        case GLUT_KEY_LEFT:
            setCoordinates(scale, theta, dx - 0.1 * scale * cosTheta, dy + 0.1 * scale * sinTheta);
            break;
        case GLUT_KEY_RIGHT:
            setCoordinates(scale, theta, dx + 0.1 * scale * cosTheta, dy - 0.1 * scale * sinTheta);
            break;
        case GLUT_KEY_DOWN:
            setCoordinates(scale, theta, dx - 0.1 * scale * sinTheta, dy - 0.1 * scale * cosTheta);
            break;
        case GLUT_KEY_UP:
            setCoordinates(scale, theta, dx + 0.1 * scale * sinTheta, dy + 0.1 * scale * cosTheta);
            break;
        default:
            break;
    }
}

void mouseFunc(int button, int state, int x, int y) {
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
        // fprintf(stderr, "\nClicked at (%d, %d)\n", x, y);
        if (transform) {
            mouse_x_down = x;
            mouse_y_down = y;
            mouse_state = state;
        }
        else {
            viewX = viewX + (2 * x / (double)windowW - 1) / viewScale;
            viewY = viewY - (2 * y / (double)windowH - 1) / viewScale;
        }
	}

	if (button == GLUT_LEFT_BUTTON && state == GLUT_UP) {
        mouse_state = state;
    }
}

void motionFunc(int x, int y) {
    mouse_x_up = x;
    mouse_y_up = y;
}

void reshape(int w, int h)
{
    windowW = w;
    windowH = h;

    windowW2 = windowW / 2.;
    windowH2 = windowH / 2.;
    invW2 = 2. / windowW;
    invH2 = 2. / windowH;

    fprintf(stderr, "\nNew size = (%d, %d)\n", w, h);
}

int main(int argc, char **argv) {
    prepare();
    readSamples();
    makeColourmap();
    particleX = size_x / (double)particleCountSqrt + 0.99;
    particleY = size_y / (double)particleCountSqrt + 0.99;
    initParticles();
    // readData();
    
	glutInit( &argc, argv );
    glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE );
    glutInitWindowSize( windowW, windowH );
    glutCreateWindow( "Hello World" );
    glutDisplayFunc( display );
    
    glutDisplayFunc(&display);
    glutIdleFunc(&display);
    glutKeyboardFunc(&keyPressed);
    glutSpecialFunc(&specialKeyPressed);
    glutMouseFunc(mouseFunc);
    glutMotionFunc(&motionFunc);
    glutReshapeFunc(&reshape);
    
    step();
    display();
    glutMainLoop();

    return 0;
}