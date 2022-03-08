/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "../Headers/OLH.h"
#include <cstring>
#include <memory>

InputOLH::InputOLH(float e_eps, uint domain) : domain(domain)
{
    g = (domain == 2) ? 2 : (floor(e_eps + 1.0));

    flip_prob = ((e_eps) / (e_eps + g - 1.0));
    bern_dist = std::bernoulli_distribution(flip_prob); // This is a prob. of generating false.
    uni_dist = std::uniform_int_distribution<uint32_t>(0, g - 1);
    pert_olh = std::vector<float>(domain, 0.0);
    keystrcache = std::vector<char *>(domain);
    char *cptr = NULL;
    lencache = new uint64_t[domain];
    std::string str = "";
    dom_index = 0;
    while (dom_index < domain)
    {
        str = std::to_string(dom_index);
        //cout << str << endl ;
        cptr = new char[str.size() + 1];
        std::strcpy(cptr, str.c_str());
        keystrcache[dom_index] = cptr;
        lencache[dom_index] = strlen(cptr);
        dom_index++;
    }
}

InputOLH::~InputOLH()
{
    for (uint i = 0; i < keystrcache.size(); i++)
        delete[] keystrcache[i];
    keystrcache.erase(keystrcache.begin(), keystrcache.end());
    delete[] lencache;
}

void InputOLH::correction()
{
    dom_index = 0;
    while (dom_index < domain)
    {
        pert_olh[dom_index] = (pert_olh[dom_index] - (1.0 * population / g)) / (flip_prob - (1.0 / g));
        dom_index++;
    }
    divBySum<vector<float>>(pert_olh);

    //cout << estimate << endl;
}

void InputOLH::perturb(const uint &key)
{
    key_str = keystrcache[key];
    hash1 = XXHash32::hash(key_str, lencache[key], population) % g;
    hash = hash1;
    if (!bern_dist(generator))
    {
        do
        {
            hash = uni_dist(generator);
        } while (hash == hash1);
    }
    dom_index = 0;
    while (dom_index < domain)
    {
        key_str = keystrcache[dom_index];
        if (hash == (XXHash32::hash(key_str, lencache[dom_index], population) % g))
            pert_olh[dom_index] += 1.0f;
        dom_index++;
    }
    population += 1;
}
