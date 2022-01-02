#include <cstdlib>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>
#include <limits.h>
#include <thread>
#include <iostream>
#include <unistd.h>

#include <cstdio>
#include <math.h>

#define GL_SILENCE_DEPRECATION
#include <GL/glew.h>
#include <GLUT/glut.h>

using namespace std;

FILE *fopen(), *f;

static unsigned int color[3];
double colorlist[999][4];

int kmax = 800;
int escape = 5;

const unsigned int W = 1920;
const unsigned int H = 1080;

GLubyte* PixelBuffer = new GLubyte[W * H * 3];

double xc = -0.5;
double yc = 0;
double step = 3.;
double zoom = 2;
double walk = 0.1;

int converge(double a, double b) {
	int k;
	double a0 = a;
	double b0 = b;
	double ap, bp, test;
	
	for(k = 0; k < kmax; k++) {
		ap = a;
		bp = b;
		
		a = ap*ap - bp*bp + a0;
		b = 2*ap*bp + b0;
		
		test = (a*a + b*b);
		
		if (test > escape) {
			break;
		}
	}
	if (k == kmax) {
		return 0;
	}
	else {
// 		return fabs((k + 1 - log(log(sqrt((double) test)))/log((double) escape)));
		return k + 1 - log( log(sqrt(test)) / log((double) escape) ) / log(2);
	};
}

void init_colorlist() {
	int r, g, b, k;
	/*
	for(k = 0; k < 256; k++) {
		r = 0;
		g = k;
		b = 255 - k;
		colorlist[k][0] = r;
		colorlist[k][1] = g;
		colorlist[k][2] = b;
	}
	for(k = 0; k < 256; k++) {
		r = k;
		g = 255 - k;
		b = 0;
		colorlist[k+256][0] = r;
		colorlist[k+256][1] = g;
		colorlist[k+256][2] = b;
	}
	for(k = 0; k < 256; k++) {
		r = 255 - k;
		g = 0;
		b = k;
		colorlist[k+512][0] = r;
		colorlist[k+512][1] = g;
		colorlist[k+512][2] = b;
	}
	*/
	
	r = 0;
	g = 0;
	b = 64;
	for(k = 0; k < 128; k++) {	
		colorlist[k][0] = r;
		colorlist[k][1] = g;
		colorlist[k][2] = b;
		
		b += 1;
		g += 2;
	}
	
	g -=1;
	//fprintf(stderr, "r: %d, g: %d, b: %d\n", r, g, b);
	for(k = 128; k < 256; k++) {		
		colorlist[k][0] = r;
		colorlist[k][1] = g;
		colorlist[k][2] = b;
		
		b -= 1;
		if (k%2 == 0) {
			b -= 1;
			g -= 1;
			r += 1;
		}
	}
	//fprintf(stderr, "r: %d, g: %d, b: %d\n", r, g, b);
	for(k = 256; k < 384; k++) {		
		colorlist[k][0] = r;
		colorlist[k][1] = g;
		colorlist[k][2] = b;
		
		g -= 1;
		if (k%2 == 0) {
			b += 1;
			g -= 1;
			r-=1;
		}
	}
	//fprintf(stderr, "r: %d, g: %d, b: %d\n", r, g, b);
}

void mandelquad(int quad) {
	int i, j, pos, count;
	
	int quadX = quad / 2;
	int quadY = quad % 2;
	
	double x;
	double step_x = step;
	double delta = step_x / W;
	double step_y = delta * H;
	
	double x0 = xc - step_x / 2. * (1-quadX);
	double y = yc - step_y / 2. * (1-quadY);
	
	int jmin = quadY * (H/2);
	int jmax = H/2 + quadY * (H/2);
	
	int imin = quadX * (W/2);
	int imax = W/2 + quadX * (W/2);
	
	for(j=jmin; j<jmax; j++) {
		x = x0;
		for(i=imin; i<imax; i++) {
			count = converge(x, y);
			//fprintf(stderr, "%.2f,%.2f,%d;",x,y,count);
			
			if (count == 0) {
				color[0] = 0;
				color[1] = 0;
				color[2] = 0;
			}
			else {
				count = count % 382+1;
				color[0] = colorlist[count][0];
				color[1] = colorlist[count][1];
				color[2] = colorlist[count][2];
			}
			
			pos = (i + j*W) * 3;
			
			PixelBuffer[pos] = color[0];
			PixelBuffer[pos+1] = color[1];
			PixelBuffer[pos+2] = color[2];
			
			x += delta;
		}
		
		y += delta;
	}
	
	//fprintf(stderr, "quad %d\n", quad);
	
	//pthread_exit(NULL);
}

void mandelbrot() {
	int i;
	//pthread_t threads[4];
	std::thread *tt = new std::thread[3];
	
	fprintf(stderr, "Render start...");
	time_t start, end;
	time(&start);
	
	for (i=0; i<3; i++) {
		tt[i] = std::thread(mandelquad, i);
	}
	
	mandelquad(3);
	
	for (int i = 0; i < 3; ++i) {
		tt[i].join();
	}
	
	//pthread_exit(NULL);
	
	PixelBuffer[((W/2 + H/2*W) * 3)] = 255;
	PixelBuffer[((W/2 + H/2*W) * 3)+1] = 0;
	PixelBuffer[((W/2 + H/2*W) * 3)+2] = 0;
	
	time(&end);
	//std::cout << difftime(end, start) << " seconds" << std::endl;
	double diff_t = difftime(end, start);
	fprintf(stderr, "Execution time = %f\n", diff_t);
	
	fprintf(stderr, "Render done\n");
	
	delete [] tt;
	
	return;
}

void key_pressed(unsigned char key, int x, int y) {
    switch (key)
    {
        case 'w':
        	step /= zoom;
        	break;
        case 's':
        	step *= zoom;
        	break;
        case 'a':
        	escape /= 2;
        	break;
        case 'd':
        	escape *= 2;
        	break;
        case 'z':
        	kmax /= 2;
        	break;
        case 'x':
        	kmax *= 2;
        	break;
        case 'r':
			fprintf(stderr, "xc = %.15f; yc = %.15f; step = %.15f; escape = %d; kmax = %d;\n", xc, yc, step, escape, kmax);
        	mandelbrot();
        	break;
        case 'q':
        	exit(0);
            break;
        // Add any keys you want to use, either for debugging or gameplay.
        default:
        	break;
    }
}

void SpecialKeys(int key, int x, int y) {
	switch (key)
	{
		case GLUT_KEY_LEFT:
			xc -= walk*step;
			break;
		case GLUT_KEY_RIGHT:
			xc += walk*step;
			break;
		case GLUT_KEY_UP:
			yc += walk*step;
			break;
		case GLUT_KEY_DOWN:
			yc -= walk*step;
			break;
	}
}

void display() {
    glClearColor( 0, 0, 0, 1 );
    glClear( GL_COLOR_BUFFER_BIT );
    glDrawPixels( W, H, GL_RGB, GL_UNSIGNED_BYTE, PixelBuffer );
    glutSwapBuffers();
}

int main( int argc, char **argv ) {
	/*xc = 0.374658203125000; yc = 0.095654296875000; step = 0.001464843750000; escape = 5; kmax = 800;
	xc = 0.374836730957031; yc = 0.095646286010742; step = 0.000011444091797; escape = 5; kmax = 800;
	xc = 0.374839019775390; yc = 0.095643404987793; step = 0.000005722045899*1.4; escape = 20; kmax = 25600;
	xc = 0.374839144945144; yc = 0.095639694956283; step = 0.000000125169754; escape = 20; kmax = 6400;*/
// 	xc = 0.374840221405029; yc = -0.095643324879150; step = 0.000008010864259; escape = 20; kmax = 25600;
// 	xc = -1.781207543141058; yc = 0.000006049649602; step = 0.000003186249265; escape = 10; kmax = 25600;
	xc = -1.401084040477871; yc = 0.000023315660655; step = 0.000000001396984; escape = 10; kmax = 51200;
	init_colorlist();
	mandelbrot();
	
	glutInit( &argc, argv );
    glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE );
    glutInitWindowSize( W, H );
    glutCreateWindow( "Mandelbrot" );
    glutDisplayFunc( display );
    
    // Bind all GLUT events to callback function.
    glutDisplayFunc(&display);
    glutIdleFunc(&display);
    glutKeyboardFunc(&key_pressed);
    glutSpecialFunc(SpecialKeys);
    
    glutMainLoop();
    return 0;
}








