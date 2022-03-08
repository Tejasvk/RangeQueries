
#pragma once

#ifndef HELPER_HEADER
#define HELPER_HEADER

#include <fstream>
#include <math.h>
#include <chrono>
#include <vector>
#include <iostream> 

typedef unsigned int uint;
typedef unsigned long long ull;
typedef unsigned long ul;


//using namespace std;
using std::ostream;
using std::vector;


const float e = exp(1);

class AllMechanisms {
public:

    double TreeRR = 0.0;
    double InputRR = 0.0;
    double TrueRange = 0.0;
    double InputOLH = 0.0;
    double TreeOLH = 0.0;
    double TreeHT = 0.0;
    double TreeHTPS = 0.0;
    double TreeHaar = 0.0;
    double InputHT = 0.0;
    double InputHTPS = 0.0;
    double cnt = 0.0;
    double length = 0.0;

    AllMechanisms() {
    };

    void clear() {
        TreeRR = 0.0;
        InputRR = 0.0;
        TrueRange = 0.0;
        InputOLH = 0.0;
        TreeOLH = 0.0;
        TreeHT = 0.0;
        TreeHTPS = 0.0;
        TreeHaar = 0.0;
        InputHT = 0.0;
        InputHTPS = 0.0;
        cnt = 0.0;
        length = 0.0;
    }


};

class AllMechanisms_min {
public:

    double InputRR = 0.0;
    double HRR = 0.0;
    
    double HaarHRR = 0.0;
    
    double TreeRR = 0.0;

    double TreeHT = 0.0;
    double TreeOLH = 0.0;
    double cnt = 0.0;
    double length = 0.0;

    AllMechanisms_min() {
    };

    void clear() {
        InputRR = 0.0;
        HaarHRR = 0.0;
        HRR = 0.0;
    
        TreeRR = 0.0;
        TreeHT = 0.0;
        TreeOLH = 0.0;

        cnt = 0.0;
        length = 0.0;
    }


};
/*
#pragma pack (1)
struct MechAnswers {
    float sm_freq_ps_u = 0.0;
    //float sm_ps_nu = 0.0;
    float sm_freq_irr = 0.0;
    //float sm_rr_nu = 0.0;
    float level_olh = 0.0;
    float level_ht = 0.0;
    float level_ht_ps = 0.0;
    float level_haar_ps = 0.0;
    uint length = 0;

    void clear() {
        sm_freq_ps_u = 0.0;
        //  sm_ps_nu = 0.0;
        sm_freq_irr = 0.0;
        //sm_rr_nu = 0.0;
        level_olh = 0.0;
        level_ht = 0.0;
        level_ht_ps = 0.0;
        level_haar_ps = 0.0;
        length = 0.0;
    }

};
*/
template<class T>
float L1(T vec1, T vec2) {
    if (vec1.size() != vec2.size())
        return -1;
    float l1 = 0.0;
    uint sz = vec1.size();

    for (uint i = 0; i < sz; i++)
        l1 += abs(vec1[i] - vec2[i]);

    return l1;
}


float power(float x, float y);
float power(float y);
void rotate(double& a, double& b);
ul ilog2(ul x);



vector<float> getWeights(int h);



template<typename T> ostream& operator<<(ostream& out, const vector<T>& v) {
    out  <<"[ ";
    size_t last = v.size() - 1;
    for (size_t i = 0; i < v.size(); ++i) {
        out << v[i];

        if (i != last)
            out << " , ";
    }
    out << " ]";
    out << std::endl;
    return out;
}

template <class T> T divBySum(T & vec) {
    float sm = 0.0;
    //cout << "=================================" << endl;
    //cout << vec;
    auto itr = vec.begin();
    for (; itr != vec.end(); ++itr)
        sm += (*itr);

    if (sm == 0.0 || sm == 1.0)
        return vec;

    itr = vec.begin();

    for (; itr != vec.end(); ++itr)
        *itr = (*itr / sm);

    //vec2[i] =(1.0*vec[i])/ sm;
    //cout << "---------------------------" << endl;
    //cout << vec;
    return vec;

}

uint inline fast_atoi(const std::string &str, const uint &len) {
    uint val = 0;
    uint cnt = len - 1;
    for (char c : str) {
        if (c == '1')
            val += (1 << cnt);
        cnt--;
    }
    return val;
}


inline uint count_1(uint& num) {
    uint cnt = 0;
    while (num != 0) {
        num = num & (num - 1);
        cnt++;
    }
    return cnt;
}

#endif