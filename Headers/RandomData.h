#pragma once
#ifndef Random_Data_h
#define Random_Data_h
#include <random>
#include "helpers.h"

class Random_Data {
public:
    std::vector<double> weights;

    int D=4;
    
    Random_Data(size_t population, uint D,float prob=0.7);
    Random_Data(std::vector<double> weights, size_t population);  
    

    inline uint sampledata() {
        return dbn(generator);
    };

private:
    std::default_random_engine generator;
    std::discrete_distribution<uint> dbn;
    void generate_data(const float& prob);
    void generate_data2(const float& prob);
    std::ofstream treefile;

        
    size_t population;
};

#endif

