#pragma once    
#ifndef IRR_HEADER
#define IRR_HEADER
#include "helpers.h"

#include <random>
#include <vector>


class InputRR {
public:
    const float e_eps;
    //const float eps;
    uint D = 8;
    float prob = 0.0;

    
    InputRR(float e_eps, int d);
    void perturb(const uint& x);
    void perturb2(uint& x);
    void reset_object();
    
    void correction();
    void correction2();
    void correction3();
    
    std::vector<double> pert_irr;
    std::vector<double> true_input;


private:
    float item = 0.0;
    double population = 0.0;
    uint i = 0;
    std::default_random_engine generator;
    std::bernoulli_distribution bernrr;

};

inline void InputRR::perturb(const uint& x) {
    true_input[x] += 1.0;
    population += 1.0;
    //std::cout << true_input << std::endl;
}
#endif




