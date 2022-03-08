/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "../Headers/LevelHT2.h"

InputHT2::InputHT2(float e_eps, uint domain, uint g) : g(g) {
    this->domain = domain;
    if ((domain & (domain - 1)) != 0)
        //this->domain = power(2,ceil(log(domain)/ log(2.0)));
        throw std::runtime_error("domain size is not a power of 2");

    g_2 = pow(2, g);
    prob_ps = (e_eps) / (g_2 + e_eps - 1.0);
    std::cout << prob_ps << " , " << g_2 << std::endl;

    rclist = vector<uint>(g, 0);
    bern_dist_ps = std::bernoulli_distribution(prob_ps); // This is a prob. of generating false.

    uni_dist = std::uniform_int_distribution<uint>(1, this->domain - 1);
    uni_dist_g = std::uniform_int_distribution<uint>(0, g_2 - 1);


    coef_dist_mulitcoef = std::vector<double>(this->domain, 0.0);
    pert_coef_dist_multicoef = std::vector<double>(this->domain, 0.0);
    pert_iht = std::vector<double>(this->domain, 0.0);

    generator.seed(time(NULL));

    uint count = 1;

    number_of_bits = vector<uint>(domain, 0);
    //number_of_bits = new uint[this->domain];
    number_of_bits[0] = 0;
    for (uint i = 1; i< this->domain; i++) {        
        count = i;
        number_of_bits[i] = count_1(count);
    }

    uint n = 0;

    coef_comb.reserve(2 * g);

    for (uint i = 0; i < g_2; i++) {
        coef_comb.erase(coef_comb.begin(), coef_comb.end());

        count = 0;
        n = i;
        vector<float> comb;
        while (n != 0) {
            coef_comb.append(n % 2 == 0 ? "1" : "-1");
            comb.push_back((n % 2 == 0 ? 1.0f : -1.0f));
            n /= 2;
            count++;
        }

        while ((g - count) > 0) {
            coef_comb.append("1");
            //comb.insert(comb.begin(), 1.0);
            comb.push_back(1.0f);
            count++;
        }
        coef_comb_float.push_back(comb);

        //cout << i << " , " << coef_comb  << endl;
        bin_comb_map[coef_comb] = i;

    }
    /*
    for (auto it : bin_comb_map){
        cout << "=========" << endl;
        std::cout << " " << it.first  << " : " << it.second << endl;
    cout << coef_comb_float[it.second] ;
    }
     * */
}

void InputHT2::perturb(const uint& key) {
    //cout << "=================" << endl;
    //    coef_comb.clear();
    coef_comb.erase(coef_comb.begin(), coef_comb.end());

    uint i = 0;
    for (; i < g; i++) {
        rc = uni_dist(generator);
        rclist[i] = rc;
        coef_dist_mulitcoef[rc] += 1.0;
        bitsagree = key & rc;
        bitsagree = number_of_bits[bitsagree];
        coef = ((bitsagree & 1) == 0) ? "1" : "-1";
        coef_comb.append(coef);
    }
    order = bin_comb_map[coef_comb];
    orderp = order;

    //if (bin_comb_map.find(coef_comb) == bin_comb_map.end())
    //    cout << coef_comb << endl;
    // cout << rclist;
    // cout << g << " , " << g_2  <<"," <<coef_comb << " , " << bin_comb_map[coef_comb] << " , " << coef_comb_float[order]; 
    if (!bern_dist_ps(generator)) {
        do {
            orderp = uni_dist_g(generator);
        } while (orderp == order);

    }
    for (i = 0; i < g; i++) {
        pert_coef_dist_multicoef [rclist[i]] += coef_comb_float[orderp][i];
    }

    // cout << coef_comb << " , " << order << " , " << coef_comb_float[order] << endl;

}

void InputHT2::correction() {
    //cout << "=========================" << endl;
    //float corr = (prob_ps * g_2 - 1.0) / (g_2 - 1.0);
    for (uint i = 1; i < domain; i++) {
        pert_coef_dist_multicoef[i] /= (prob_ps * std::max(1.0, coef_dist_mulitcoef[i]));
    }
    pert_coef_dist_multicoef[0] = 1.0;

    fwht();

    double sm_iht_mult = 0.0;
    uint i = 0;
    for (; i < domain; i++) {
        sm_iht_mult += pert_coef_dist_multicoef[i];
    }
    for (i = 0; i < domain; i++) {
        pert_iht[i] = pert_coef_dist_multicoef[i] / sm_iht_mult;
    }

    //   std::cout << pert_iht_multi_coef;

}

/**
 * Fast Walsh-Hadamard transform
 */

void InputHT2::fwht(void) {
    const ul l2 = ilog2(pert_coef_dist_multicoef.size()) - 1;

    for (ul i = 0; i < l2; ++i) {
        for (ul j = 0; j < (1 << l2); j += 1 << (i + 1)) {
            for (ul k = 0; k < (1 << i); ++k) {
                //rotate(data[j + k], data[j + k + (1 << i)]);
                rotate(pert_coef_dist_multicoef[j + k], pert_coef_dist_multicoef[j + k + (1 << i)]);
                //cout << "hii" << endl;
            }
        }
        //      cout << i << endl;
    }
    //std::cout << "inverted "  << data;
}
