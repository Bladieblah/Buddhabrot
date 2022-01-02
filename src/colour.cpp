#include <vector>
#include <iostream>

#include "colour.hpp"
#include "interp.hpp"

using namespace std;

Colour::Colour(vector<float> _x, vector< vector<float> > _y, int _size) {
    int i;
    
    size = _size;
    
    Interp1d interp(_x, _y);
    
    float p = 0;
    float dp = 1. / (size - 1.);
    
    for (i=0; i<size; i++) {
        vector<float> result = interp.getValue(p);
        map.push_back(result);
        
        p += dp;
    }
}

vector<float> Colour::get(float p) {
    int i = (int)(p * size);
    
    return map[i];
}

void Colour::apply(float *colourMap) {
    for (int i=0; i<size; i++) {
        for (int j=0; j<3; j++) {
            colourMap[3 * i + j] = map[i][j];
        }
    }
}

void Colour::apply(uint32_t *colourMap) {
    for (int i=0; i<size; i++) {
        for (int j=0; j<3; j++) {
            colourMap[3 * i + j] = map[i][j] * uint32_max;
        }
    }
}