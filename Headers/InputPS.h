//#ifndef INPUTPS_H
//#define INPUTPS_H
//#include "InputPS.h"
#pragma once
#include "helpers.h"

class InputPS {
public:
    int D;
    float e_eps = 0.0;
    float prob = 0.0;
    vector <bool>bern;

    //Eigen::Matrix <float, Dynamic, Dynamic>mat;
    void perturb(int x, int p);
    int perturb(int &x);
    vector <float> correction();

    InputPS(int D, float e_eps, float population) {
        this->D = D;
        this->e_eps = e_eps;
        this->population = population;
        prob = (e_eps - 1.0) / (D + e_eps - 1.0);
        bern_dist = std::bernoulli_distribution(prob); // This is a prob. of generating false.
        uni_dist = std::uniform_int_distribution<int>(0, D);

        pert.resize(D);
        input.resize(D);
        std::fill(pert.begin(), pert.end(), 0.0);
        std::fill(input.begin(), input.end(), 0.0);        
        //cout << pert.size() << endl;
        //populate_matrix();
        generate_randomness();
    }
    ~InputPS();
private:
    float population;
    size_t freq = 0;
    int choice = 0;
    void populate_matrix();
    void generate_randomness();
    void generate_options();
    std::default_random_engine generator;
    std::bernoulli_distribution bern_dist;
    std::vector<float>options;
    std::vector<float> pert;
    std::vector<float> input;


    std::uniform_int_distribution<int> uni_dist;
};

InputPS::~InputPS() {
    //pert.erase(pert.begin(),pert.end());
    //input.erase(input.begin(),input.end());
}

vector<float> InputPS::correction() {
    vector<float> pert_e = divBySum<vector<float>>(pert);
    vector<float> true_e = divBySum<vector<float>>(input);
    //cout <<endl<< "correcting ps" << endl;
    vector<float> corrected(D);
    for (auto i = 0; i < D; i++)
        corrected[i] = ((pert_e[i] * D) + prob - 1.0) / (prob * (D + 1.0) - 1.0);
    corrected = divBySum<vector<float>>(corrected);
    
    cout << true_e << endl;
    cout << corrected << endl;
    cout << "L1=" << L1(corrected,true_e) << endl;
    return corrected;
}

void InputPS::generate_randomness() {
    bool flag = false;
    uni_dist.reset();
    for (auto p = 0; p < population; p++) {
        flag = bern_dist(generator);
        bern.push_back(flag);
        choice = uni_dist(generator);
        options.push_back(choice);
    }
}

void InputPS::perturb(int x, int p) {
    input[x] += 1.0;
    if (bern[p]) {
        pert[x] += 1.0;
    } else {
        pert[options[p]] += 1.0;
    }
}

/*void InputPS::populate_matrix() {
    mat.resize(D, D);
    for (auto i = 0; i < D; i++) {
        for (auto j = 0; j < D; j++) {
            if (i == j)
                mat(i, j) = (prob + 1.0 / D);
            else
                mat(i, j) = (1.0 - prob) / (D);
        }
    }
}*/