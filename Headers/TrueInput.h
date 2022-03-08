#pragma once
#ifndef TRUEINPUT_H
#define TRUEINPUT_H
#include <vector>

class TrueInput {
public:
    const uint D;
    std::vector<double> true_input;
    void reset_object();
    TrueInput(uint d) : D(d) {
        true_input = vector<double>(d, 0.0);
    }

    void normalize();
    

private:
};


#endif