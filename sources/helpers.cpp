/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "../Headers/helpers.h"


float power(float x, float y) {
    return pow(x, y);
}

float power(float y) {
    return pow(e, y);
}

vector<float> getWeights(int h) {
    vector<float> weights;
    int height = pow(2, h);
    for (int i = 0; i < height; i++)
        weights.push_back(rand() % height);
    //cout << weights;
    return weights;
}


void rotate(double& a, double& b) {
    static double t;	
    t = a;
    a = a + b;
    b = t - b;
}

// Integer log2

ul ilog2(ul x) {
    ul l2 = 0;
    for (; x; x >>= 1) ++l2;
    return l2;
}


