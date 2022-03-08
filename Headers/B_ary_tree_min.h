/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   B_ary_tree_min.h
 * Author: godzilla
 *
 * Created on 28 March 2018, 12:06
 */

#pragma once
#ifndef B_ARY_TREE_MIN_HEADER
#define B_ARY_TREE_MIN_HEADER
#include <tuple>
#include <queue>
#include <stack>
#include <memory>
#include "node_b_l.h"
#include "InputRR.h"
#include "LevelHT.h"
#define  OLH 0
#if OLH == 1
#include "OLH.h"
#endif

//using namespace std;

class B_ary_tree_min
{

  public:
    B_ary_tree_min(int l, int r, float e_eps, int b);

    ~B_ary_tree_min()
    {
        for (auto itr = arraytree.begin(); itr != arraytree.end(); itr++)
        {
            if (*itr)
                delete *itr;
        }

        arraytree.shrink_to_fit();
        /*
        std::vector<node_b_l *> temptree;
    
        arraytree.swap(temptree);

        arraytree.erase(arraytree.begin(), arraytree.end());
        arraytree.shrink_to_fit();*/
        std::cout << "done inside tree destructor" << std::endl;
    }
    std::vector<node_b_l *> arraytree;
    AllMechanisms_min answer;

    void buildArrayTree();
    // std::vector<int> getPath(const uint leaf);
    void correction();
    void constrained_Inference();
    void perturb(const uint &x);
    void reset_object();
    void eval_range_itr(vector<int> &range);

    AllMechanisms_min evalRange(vector<int> &range);

  private:
    const int l_range;
    const int r_range;
    const float e_eps;
    const int b;

    int level = 1, height;
    uint trueindex = 0;
    uint noisy = 0;
    uint i = 0;
    uint j = 0;
    uint leafoffset = 0;

    float prob_rr = 0.0;
    float item = 0.0;

    void weighted_averaging(node_b_l *root);
    void mean_Consistency(node_b_l *root);

    void weighted_averaging_itr(node_b_l *root);
    void mean_Consistency_itr(node_b_l *root);

    void buildIndexCache();

    node_b_l *n;
    AllMechanisms_min c_sum;
#if OLH == 1
    std::vector<std::unique_ptr<InputOLH>> olharray;    
#endif
    std::vector<std::unique_ptr<InputHT>> lhtarray;
    std::vector<std::unique_ptr<InputRR>> uirrarray;

    std::stack<std::tuple<std::vector<int>, node_b_l *>> nodeStack;
    std::vector<std::vector<uint>> trueindexcache;

    std::vector<float> phi = {0.0, 0.0};
    vector<int> phi_range = {0, 0};
    vector<int> phi2 = {0, 0};

    std::default_random_engine generator;
    std::uniform_int_distribution<int> uni_level_sampler; // uni_level_sampler_dht;
};

#endif