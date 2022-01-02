#include <vector>

#include "interp.hpp"

using namespace std;

Interp1d::Interp1d(vector<float> _x, vector< vector<float> > _y) {
    for (int i=0; i<(int)_x.size(); i++) {
        x.push_back(_x[i]);
        
        vector<float> temp;
        for (int j=0; j<(int)_y[i].size(); j++)
            temp.push_back((float)_y[i][j] / 255.);
        
        y.push_back(temp);
    }
}

float Interp1d::interpolate(float x0, float x1, float y0, float y1, float p) {
    return y0 + (y1-y0) / (x1-x0) * (p-x0);
}

vector<float> Interp1d::getY(int i) {
    int j;
    vector<float> result;
    
    for (j=0; j<(int)y[i].size(); j++)
        result.push_back(y[i][j]);
    
    return result;
}

vector<float> Interp1d::getValue(float p) {
    int i, j;
    vector<float> result;
    
    if (p >= x[(int)x.size() - 1])
        return getY((int)x.size() - 1);
    
    for (i=0; i<(int)x.size(); i++) {
        if (p < x[i])
            break;
    }
    
    if (i == 0)
        return getY(i);
    
    for (j=0; j<y[i].size(); j++) {
        result.push_back(interpolate(x[i-1], x[i], y[i-1][j], y[i][j], p));
    }
    
    return result;
}