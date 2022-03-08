
#pragma once
#ifndef LEVEL_HT_HEADER
#define LEVEL_HT_HEADER
#include "helpers.h"
#include <iostream>

#include <random>


class InputHT {
public:
    //~InputHT();

    InputHT(float e_eps, uint domain,bool flag=true);
    void perturb(const int& key,bool sign=true);
    void correction();
    void reset_object();
    std::vector<double> pert_iht;
    std::vector<double>pert_coef_dist;
    
    std::vector<float> pert_dht_ht;
private:
    uint rc, bitsagree = 0;
    double coef, prob;
    const uint domain;
    const bool signal_sign;
    std::uniform_int_distribution<uint> uni_dist;
    std::bernoulli_distribution bern_dist ;
    std::vector<double>coef_dist;
    std::vector<uint>number_of_bits;

    void fwht();

    std::default_random_engine generator;



};

#endif
