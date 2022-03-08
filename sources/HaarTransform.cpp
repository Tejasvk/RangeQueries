/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "../Headers/HaarTransform.h"
#include <iostream>
#include <exception>
#include <ctime>
#define debug 0
#if debug == 1
#include <algorithm> // std::count_if
#endif

HaarTransform::HaarTransform(float e_eps, uint domain) : domain(domain) , height(log2(domain)), e_eps(e_eps)
{

    if ((domain & (domain - 1)) != 0)
        throw std::runtime_error("domain size is not a power of 2");

    
    coef_set.reserve(2 * height + 1);
    this->prob = e_eps / (1.0 + e_eps);
    bern = std::bernoulli_distribution(prob);

    levelfreq = vector<double>(height, 0.0);
    pert_haar_ht = vector<double>(domain, 0.0);

    haar_transform();
    first_coef=0.0;
    uni_level_sampler_dht = std::uniform_int_distribution<int>(0, height - 1);

    generator.seed(time(NULL));

    ht_array.reserve(height);
    uint count_level = 1;
    while (count_level < height)
    {
      //  std::cout << count_level << " , " << pow(2, count_level) << "\n";
        ht_array.push_back(std::make_unique<InputHT>(e_eps, (uint)pow(2, count_level), false));
        count_level++;
    }
}

HaarTransform::~HaarTransform()
{
    std::cout << "inside Haar destructor!" << std::endl;
    ht_array.clear();
    for (uint i = 0; i < domain; i++)
        haar_coef_index_list[i].shrink_to_fit();

    pert_haar_ht.shrink_to_fit();
    levelfreq.shrink_to_fit();

    std::cout << "deleted Haar object!" << std::endl;
}

/*
We reset counters for all the objects but retains the caches.
Use this function when you want to reuse the object after clearing the counters.
*/
void HaarTransform::reset_object()
{

    ht_array.clear();
    
    ht_array.erase(ht_array.begin(), ht_array.end());

    uint count_level = 1;    
    while (count_level < height)
    {
        ht_array.push_back(std::make_unique<InputHT>(e_eps, (uint)pow(2, count_level), false));
        count_level++;
    }
    std::fill(levelfreq.begin(), levelfreq.end(), 0.0f);
    std::fill(pert_haar_ht.begin(), pert_haar_ht.end(), 0.0f);

    population = 0.0;
    first_coef = 0.0f;
    std::cout << "haar tree reset" << std::endl;
}


void HaarTransform::haar_transform()
{
    // We first see if Haar transform matrix is stored locally. If so, we load it in the memory. Otherwise, we have to compute it.
    // Note that precomputing (signed) signal indices for each item, level pair is a very compute intense operation. 
    std::string treefilename = "../treefiles/haar_tree_" + std::to_string(domain) + ".csv";
    std::ifstream in_treefile(treefilename);
    
    if (in_treefile.good())
    {
        
        std::string line;
        if (in_treefile.is_open())
        {
            clock_t begin = clock();
            std::string token;
            while (getline(in_treefile, line))
            {
                //std:: cout << "=============== \n";
                //std::cout << line << "\n";
                std::vector<int> list;
                list.reserve(height);
                token.reserve(5);

                for (char c : line)
                {
                    if (c != ',')
                        token += c;
                    else
                    {
                        list.emplace_back(stoi(token));
                        token.erase(token.begin(), token.end());
                    }
                }
                haar_coef_index_list.emplace_back(list);
            }
            in_treefile.close();
            clock_t end = clock();
            std::cout << "haar tree loaded in " << double(end - begin) / CLOCKS_PER_SEC << " seconds \n";
        }
    }
    else
    {
        // We compute and store the Haar transform matrix for each item. 
        uint n = domain;
        vector<float> avg = vector<float>(domain, 0.0f);
        vector<float> diff = vector<float>(domain, 0.0f);
        vector<float> xc = vector<float>(domain, 0.0f);
        //float xc[domain];
        uint count_level = 1;
        uint index = 0;
        //uint d = log2(domain);
        //size_t szoffloat = sizeof(float);
        haar_coef_index_list.reserve(domain);
        for (uint dom = 0; dom < domain; dom++)
        {

            std::fill(xc.begin(), xc.end(), 0.0f);
            //memset(xc, 0.0f, szoffloat * domain);

            n = domain;
            if (dom % 10000 == 0)
                std::cout << dom << std::endl;

            xc[dom] = 1.0f;
            while (n > 1)
            {
                for (uint i = 0; i < (n >> 1); i++)
                {
                    avg[i] = (xc[2 * i] + xc[2 * i + 1]) / 2.0f;
                    diff[i] = xc[2 * i] - avg[i];
                }
                for (uint i = 0; i < (n >> 1); i++)
                {
                    xc[i] = avg[i];
                    xc[i + (n >> 1)] = diff[i];
                }
                n = n >> 1;
            }
           // std::cout << xc;
            count_level = 1;
            std::vector<int> list;
            list.reserve(height);
            while (count_level <= height)
            {
                leftend = 1 << (count_level - 1);
                rightend = 1 << count_level;
                index = 0;
                while (leftend < rightend)
                {
                    if (xc[leftend] > 0.0f)
                        list.emplace_back(index + 1);
                    else if (xc[leftend] < 0.0f)
                        list.emplace_back(-(index + 1));
                    index++;
                    leftend++;
                }
                count_level += 1;
            }
            haar_coef_index_list.emplace_back(list);

        }
        clock_t begin = clock();

        std::ofstream treefile;
        treefile.open(treefilename);
        for (uint dom = 0; dom < domain; dom++)
        {
            index = 0;
            for (auto d : haar_coef_index_list[dom])
            {
                treefile << d;
                if (index != domain)
                    treefile << ",";
                index++;
            }
            treefile << "\n";
        }

        treefile.close();
        clock_t end = clock();
        std::cout << "haar tree saved in " << double(end - begin) / CLOCKS_PER_SEC << " seconds \n";

    }
    std::cout << "haar  tree allocated" << std::endl;
}

/*
Sample a level and take Hadamard transform and then perturb a sampled coefficient.
*/
void HaarTransform::perturb(const uint &x)
{

    level = uni_level_sampler_dht(generator);
    population += 1.0;
    levelfreq[level] += 1.0;
 
    signal = haar_coef_index_list[x][level];
    sign = (signal > 0);

    signal = (signal > 0) ? signal - 1 : signal + 1;
    if (level > 0){
        ht_array[level - 1]->perturb(signal, sign);
    }else {
         coef = (sign) ? 1.0 : -1.0;
         if (!bern(generator))
            coef = -coef;
        first_coef += coef;
    }
   
}


/*
Correct the haar tree levelwise.
*/
void HaarTransform::correction()
{

    uint count = 2;
    coef = 0.0;
    uint index = 0;

    while (count <= height)
    {
        leftend = 1 << (count - 1);
        rightend = 1 << count;
        coef = pow(0.5, height - count + 1);
        ht_array[count - 2]->correction();

        index = 0;
        while (leftend < rightend)
        {
            pert_haar_ht[leftend] = (population / levelfreq[count-1]) * coef * ht_array[count - 2]->pert_coef_dist[index];
            index++;
            leftend++;
        }
        count++;
    }
    // The 0th coefficient is always hardcoded to N/D. 
    pert_haar_ht[0] = (population / (double) domain);

    coef = pow(0.5, height);        
    pert_haar_ht[1] = (population/levelfreq[0])*coef*first_coef/(2.0*prob-1.0);  

}

#if debug == 1
bool Iszero(double i) { return (i != 0.0); }
#endif

/*
We iteratively compute a given range by aggregating over <= 2h + 1 coefficients.
For more details, please consult 
https://www.semanticscholar.org/paper/The-Discrete-Wavelet-Transform-and-Wavelet-Synopses-Garofalakis/0560f3446fe329dd12d07ee93332e4fef13b665c 
*/
double HaarTransform::eval_range_itr(std::vector<int> const &range)
{
    leftend = floor((range[0] + domain) / 2);
    rightend = floor((range[1] + domain) / 2);
    qval_ht = 0.0;
#if debug == 1
    //std::cout << "================" << std::endl;
    vector<doubles> weights_vec(height + 1, 0.0);
    no_nodes = 0;
    //std::cout << range;
#endif
    weights = 0.0;
    count = 1;

    while (leftend > 0)
    {

        weights = 0.0;
        subtree[0] = (1 << count) * (leftend)-domain;
        subtree[1] = subtree[0] + floor(((1 << count) - 1) / 2);
        overlap[0] = std::max(subtree[0], range[0]);
        overlap[1] = std::min(subtree[1], range[1]);
        if (overlap[0] <= overlap[1])
        {
            weights = overlap[1] - overlap[0] + 1.0;
#if debug == 1
            weights_vec[height - count + 1] = overlap[1] - overlap[0] + 1.0f;
#endif
        }

        subtree[1] = (1 << count) * (leftend + 1) - 1 - (domain);
        subtree[0] = subtree[1] - floor(((1 << count) - 1) / 2);
        overlap[0] = std::max(subtree[0], range[0]);
        overlap[1] = std::min(subtree[1], range[1]);
        if (overlap[0] <= overlap[1])
        {
            weights -= (overlap[1] - overlap[0] + 1.0);
#if debug == 1
            weights_vec[height - count + 1] -= overlap[1] - overlap[0] + 1.0f;
#endif
        }
#if debug == 1
        if (weights != 0.0f)
            no_nodes++;
#endif

        qval_ht += (pert_haar_ht[leftend] * weights);

        coef_set.insert(leftend);
        leftend = floor(leftend / 2);
        count++;
    }
    weights = range[1] - range[0] + 1.0;

#if debug == 1
    weights_vec[0] = range[1] - range[0] + 1.0f;
    int mycount = count_if(weights_vec.begin(), weights_vec.end(), Iszero);
  //  std::cout << mycount << " left end " << weights_vec;
#endif

    qval_ht += (pert_haar_ht[leftend] * weights);
    count = 1;
#if debug == 1
    std::fill(weights_vec.begin(), weights_vec.end(), 0.0);
#endif
    while (rightend > 0)
    {

        if (coef_set.count(rightend) == 0)
        {
            weights = 0.0;
            subtree[0] = (1 << count) * (rightend)-domain;
            subtree[1] = subtree[0] + floor(((1 << count) - 1) / 2);
            overlap[0] = std::max(subtree[0], range[0]);
            overlap[1] = std::min(subtree[1], range[1]);
            if (overlap[0] <= overlap[1])
            {
#if debug == 1
                weights_vec[height - count + 1] = overlap[1] - overlap[0] + 1.0f;
#endif
                weights = overlap[1] - overlap[0] + 1.0;
            }
            subtree[1] = (1 << count) * (rightend + 1) - 1 - (domain);
            subtree[0] = subtree[1] - floor(((1 << count) - 1) / 2);
            overlap[0] = std::max(subtree[0], range[0]);
            overlap[1] = std::min(subtree[1], range[1]);
            if (overlap[0] <= overlap[1])
            {
#if debug == 1
                weights_vec[height - count + 1] -= (overlap[1] - overlap[0] + 1.0f);
#endif
                weights -= (overlap[1] - overlap[0] + 1.0);
            }
            qval_ht += (pert_haar_ht[rightend] * weights);
#if debug == 1
            if (weights != 0.0f)
                no_nodes++;
#endif

            coef_set.insert(rightend);
        }
        rightend = floor(rightend / 2);
        count++;
    }
#if debug == 1
    mycount = count_if(weights_vec.begin(), weights_vec.end(), Iszero);
  //  std::cout << mycount << " right end " << weights_vec;
#endif
    std::unordered_set<int>().swap(coef_set);
    return qval_ht / population;
}
