/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "../Headers/RandomData.h"
Random_Data::Random_Data(size_t population, uint D, float prob)
{

    generator.seed(time(NULL));

    this->D = D;
    this->population = population;
    generate_data2(prob);
    //cout << D;

    dbn = std::discrete_distribution<uint>(this->weights.begin(), this->weights.end());
    weights = std::vector<double>();
    weights.shrink_to_fit();
}

Random_Data::Random_Data(std::vector<double> weights, size_t population)
{
    this->weights = weights;
    this->population = population;
    dbn = std::discrete_distribution<uint>(this->weights.begin(), this->weights.end());
}

void Random_Data::generate_data2(const float &prob)
{
    weights = std::vector<double>(D, 0.0);
    std::cout << "Cauchy distribution prob. = " << prob << std::endl;
   // std::poisson_distribution<uint> dbn(prob);
    std::cauchy_distribution<double> dbn(D * prob, D / 10.0);
//    std::lognormal_distribution<double> distribution(0.0,1.0);
   // std::exponential_distribution<double> distribution(3.5);
    //  std::geometric_distribution<int> distribution(prob);

    //std::uniform_int_distribution<uint> dbn(0,D);
    size_t p = 0;
    uint number = 0;
    while (p <= 200000000)
    {
        number = (uint)abs(dbn(generator));
        //number = (dbn(generator));
        //std::cout << dbn(generator) << "\n";
        if ((number >= 0) && (number < D))
        {
            weights[number] += 1.0;
            p++;
        }
    }

    divBySum<vector<double>>(weights);
    //  std::cout << weights;
/*
    treefile.open("dbn_"+std::to_string(D)+".csv" ,std::ios::app);
    uint index = 0;
    for (uint dom = 0; dom < D; dom++)
    {
        treefile << weights[dom];
        if (index != D)
            treefile << ",";
        index++;
    }
    treefile << "\n";

    treefile.close();
    //  std::cout << sum << "\n";
*/
  
}

void Random_Data::generate_data(const float &prob)
{
    std::bernoulli_distribution bern(prob);
    uint cnt = 0;
    uint cnt2 = 0;
    uint x = 0;
    uint d = log2(D);
    weights = std::vector<double>(D, 0.0);

    for (size_t p = 0; p < 100000000; p++)
    {
        cnt = 0;
        x = 0;
        cnt2 = d - 1;
        while (cnt < d)
        {
            if (bern(generator))
                x = x + (1 << cnt2);
            cnt2--;
            cnt++;
        }
        weights[x] += 1.0;
    }
    //std::cout << inputdbn;

    //divBySum<std::vector<double>>(weights);
    std::cout << weights;
}
