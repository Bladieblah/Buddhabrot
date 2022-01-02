#ifndef COLOUR_H
#define COLOUR_H

#include <vector>

class Colour {
public:
    Colour(std::vector<float> _x, std::vector< std::vector<float> > _y, int _size);
    
    void apply(float *colourMap);
    void apply(uint32_t *colourMap);
    std::vector<float> get(float p);
    
    std::vector< std::vector<float> > map;
    
    int size;
private:
    uint32_t uint32_max = 4294967295;
};

#endif