#pragma once
#ifndef OLH_h
#define OLH_h   

#include "xxhash32.h"
#include "helpers.h"
#include <random>

class InputOLH {
public:
    ~InputOLH();

    InputOLH(float e_eps, uint domain);
    void perturb(const uint& key);
    void correction();
    std::vector<float> pert_olh;

private:
    uint g;
    const uint domain;
    uint population = 0;
    uint32_t hash = 0,hash1=0;
    float flip_prob = 1.0;
    char* key_str;
   
    uint64_t* lencache = NULL;
    //std::shared_ptr<uint64_t> lencache;

    uint dom_index = 0;

    std::vector<char*> keystrcache;
    std::uniform_int_distribution<uint32_t> uni_dist;
    std::bernoulli_distribution bern_dist;

    std::default_random_engine generator;

};

#endif


