
#ifndef HaarTransform_h
#define HaarTransform_h
#pragma once
#include "helpers.h"
#include <random>
#include <memory>
#include <vector>


#include "LevelHT.h"
#include <unordered_set>

class HaarTransform {
public:

    HaarTransform(float e_eps, uint domain);
    ~HaarTransform();
    
    void haar_transform();
    void correction();
    void reset_object();

    void perturb(const uint & x);
    double eval_range_itr(vector<int>const & range);
    vector<double> pert_haar_ht;
  
    
    vector<double> levelfreq;
    //vector<float> truedist;

    double qval_ht = 0.0;
    uint no_nodes= 0;

    

private:   
    bool sign = false;
    vector<std::unique_ptr<InputHT>> ht_array;
    vector<vector<int>> haar_coef_index_list;

    int subtree[2] = {0, 0};
    int overlap[2] = {0, 0};
    int leftend = 0;
    int rightend =0;
    int signal =0;
    const uint domain;
    std::unordered_set<int> coef_set;
    
    const uint height = 0;
    uint level = 1;
    uint count =1;
    
    float prob = 1.0f;
    const float e_eps;
    double weights=0.0f;
    double first_coef = 0.0;
    double coef = 1.0;
    double population = 0.0;  
    
    std::bernoulli_distribution bern;

    std::uniform_int_distribution<int> uni_level_sampler_dht;
    std::default_random_engine generator;


};
#endif






