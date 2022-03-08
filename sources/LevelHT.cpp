/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "../Headers/LevelHT.h"

InputHT::InputHT(float e_eps, uint domain, bool flag): signal_sign(flag) , domain(domain)
{
    if ((domain & (domain - 1)) != 0)
        throw std::runtime_error("domain size is not a power of 2");

    prob = e_eps / (1.0f + e_eps);
    bern_dist = std::bernoulli_distribution(prob); // This is a prob. of generating false.
    if (signal_sign)
        uni_dist = std::uniform_int_distribution<uint>(1, this->domain - 1);
    else
        uni_dist = std::uniform_int_distribution<uint>(0, this->domain - 1);

    coef_dist = std::vector<double>(this->domain, 0.0);
    pert_coef_dist = std::vector<double>(this->domain, 0.0);
    pert_iht = std::vector<double>(this->domain, 0.0);

   
    generator.seed(time(NULL));

    uint count = 1;

    number_of_bits = std::vector<uint>(domain, 0);
    number_of_bits[0] = 0;
    for (uint i = 1; i < this->domain; i++)
    {
        count = i;
        number_of_bits[i] = count_1(count);
    }
}


void InputHT::perturb(const int &key, bool sign)
{
    rc = uni_dist(generator);
    coef_dist[rc] += 1.0;
    bitsagree = abs(key) & rc;
    bitsagree = number_of_bits[bitsagree];

    coef = ((bitsagree & (uint)1) == 0) ? 1.0 : -1.0;
    if (!sign)  
        coef = -coef;

    if (!bern_dist(generator))
        coef = -coef;
    pert_coef_dist[rc] += (coef);
}
void InputHT::reset_object(){
    for (size_t i = 0; i < domain; i++)
    {
        coef_dist[i] = 0.0;
        pert_iht[i] = 0.0;
        pert_coef_dist[i] =0.0;
    }
 std::cout << "HRR object reset." << std::endl;
   
}

void InputHT::correction()
{
    double corr = (2.0 * prob - 1.0);
    for (uint i = 0; i < domain; i++)
    {
        pert_coef_dist[i] /= corr ;
        if (signal_sign)
            pert_coef_dist[i] /= std::max(1.0, coef_dist[i]);

    }
    if (signal_sign)    
        pert_coef_dist[0] = 1.0;
    
    fwht();
    if (signal_sign)
    {
        double sm_iht = 0.0;
        for (uint i = 0; i < domain; i++)
        {
            sm_iht += pert_coef_dist[i];
        }
        for (uint i = 0; i < domain; i++)        
            pert_iht[i] = pert_coef_dist[i] / sm_iht;
        
    }

}

/**
 * Fast Walsh-Hadamard transform
 */

void InputHT::fwht(void)
{
    const ul l2 = ilog2(pert_coef_dist.size()) - 1;

    for (ul i = 0; i < l2; ++i)
    {
        for (ul j = 0; j < (1 << l2); j += 1 << (i + 1))
        {
            for (ul k = 0; k < (1 << i); ++k)
            {
                rotate(pert_coef_dist[j + k], pert_coef_dist[j + k + (1 << i)]);
            }
        }
    }

}
