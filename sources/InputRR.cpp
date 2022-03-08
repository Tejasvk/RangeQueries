/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "../Headers/InputRR.h"
InputRR::InputRR(float e_eps, int d) : e_eps(e_eps), D(d)
{

    //prob = (sqrt(e_eps)) / (1.0 + sqrt(e_eps));

    prob = (1.0) / (1.0 + e_eps);

    pert_irr = std::vector<double>(D, 0.0);
    true_input = std::vector<double>(D, 0.0);
    bernrr = std::bernoulli_distribution(prob);
}

void InputRR::reset_object()
{
    population = 0.0;
    for (size_t i = 0; i < D; i++)
    {
        pert_irr[i] = 0.0;
        true_input[i] = 0.0;
    }
    std::cout << "inputRR object reset!" << std::endl;
}

/*
We simulation InputRR (OUE) using binomial distribution.
*/
void InputRR::correction()
{
 
    std::binomial_distribution<uint> bino1, bino2;
    double sum = 0.0;
    for (uint i = 0; i < D; i++)
    {
        bino1 = std::binomial_distribution<uint>(true_input[i], 0.5);
        bino2 = std::binomial_distribution<uint>(population - true_input[i], prob);
        pert_irr[i] = 0.0;
        pert_irr[i] += ((bino1(generator) + bino2(generator)));
        pert_irr[i] /= (1.0 * population);
        pert_irr[i] = (pert_irr[i] - prob) / (0.5 - prob);
    }
    
}

/*
We simulation InputRR (OUE) using binomial distribution.
Use this if you want to normalize the reconstructed counts.
*/
void InputRR::correction3()
{

    std::binomial_distribution<uint> bino1, bino2;
    double sum = 0.0;
    for (uint i = 0; i < D; i++)
    {
        bino1 = std::binomial_distribution<uint>(true_input[i], 0.5f);
        bino2 = std::binomial_distribution<uint>(population - true_input[i], prob);
        pert_irr[i] = 0.0;
        pert_irr[i] += ((bino1(generator) + bino2(generator)));
        pert_irr[i] /= (1.0 * population);
        pert_irr[i] = (pert_irr[i] - prob) / (0.5 - prob);

        sum +=pert_irr[i];
    }
    /*
    We redestribute the difference between 1 and the total sum. 
    */
    double diff = abs(1.0-sum)/D;
    for (size_t i= 0; i < D; i++)    
        pert_irr[i]+= diff;   
    
}

/*
Use this method when you can afford to perturb all D bits from each user.  
*/
void InputRR::perturb2(uint &x)
{
    i = 0;
    population += 1.0;
    while (i < D)
    {
        item = 0.0;
        if (i == x)
            item = 1.0;
        if (!bernrr(generator))
            item = 1.0 - item;
        pert_irr[i] += item;
        i++;
    }
    //   cout << true_input << endl;
}

/*
Use this when along with perturb2.
*/
void InputRR::correction2()
{
    for (i = 0; i < D; i++)
    {
        pert_irr[i] /= population;
        pert_irr[i] = (pert_irr[i] + prob - 1.0) / (2.0 * prob - 1.0);
    }
}
