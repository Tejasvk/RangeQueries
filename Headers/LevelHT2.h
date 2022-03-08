
#pragma once
#include "helpers.h"
#include <iostream>

#include <random>
#include <unordered_map>

class InputHT2 {
public:

    InputHT2(float e_eps, uint domain,uint g);  
    void perturb(const uint& key);
    void correction();
    std::vector<double> pert_iht;
    std::vector<double> pert_dht_ht;
    

private:
    uint rc, bitsagree = 0, order = 0,orderp = 0;
    const uint g;
    float prob_ps;
    std::string coef= "1",coef_comb = "";
    
    uint domain,g_2;
    std::uniform_int_distribution<uint> uni_dist;
    std::uniform_int_distribution<uint> uni_dist_g;
    std::bernoulli_distribution  bern_dist_ps;
    
    std::vector<double>coef_dist_mulitcoef;
    std::vector<double>pert_coef_dist;
    std::vector<uint>number_of_bits;
    std::unordered_map<std::string,uint> bin_comb_map;
    std::vector<uint>rclist;
 
    std::vector<vector<float>>coef_comb_float;


    std::vector<double>pert_coef_dist_multicoef;

    void fwht();
    
    std::default_random_engine generator;


};


