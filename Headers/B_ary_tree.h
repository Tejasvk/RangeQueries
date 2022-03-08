
/* 
 * File:   B_ary_tree.h
 * Author: godzilla
 *
 * Created on 28 March 2018, 12:06
 */
#pragma once
#ifndef B_ARY_TREE_HEADER
#define B_ARY_TREE_HEADER

//#include <iostream>
//#include <random>
#include <tuple>
#include <queue>
#include <stack>
#include <memory>
#include "node_b.h"
#include "InputRR.h"
#include "LevelHT.h"
#include "OLH.h"



//using std::ostream;
//using std::vector;

class B_ary_tree {
public:
    B_ary_tree(int l, int r, float e_eps, int b);

    ~B_ary_tree() {
        //for (auto itr = arraytree.begin(); itr != arraytree.end(); ++itr) {
        //    if (*itr) delete *itr;
        //}
        arraytree.erase(arraytree.begin(), arraytree.end());
        
    }
    vector<node_b *> arraytree;
    


    vector<node_b*> buildArrayTree();
    void correction();
    void constrained_Inference();
    void consistency();
    void perturb(const uint& x);
    void reset_object();
    AllMechanisms answer;
    AllMechanisms eval_range_itr(vector<int>&range);


    
private:
    const int l_range;
    const int r_range;
    const float e_eps;
    const int b;

    int level = 1;
    int height;
    int temp_phi[2]={0,0};
    uint trueindex = 0;
    uint noisy = 0;
    uint i = 0;
    uint j = 0;
    uint leafoffset = 0;
    
    float item = 0.0;
    
    void weighted_averaging(node_b *root);
    void mean_Consistency(node_b *root);
    void buildIndexCache3();   

    node_b *n;
    AllMechanisms c_sum;

    std::vector<std::unique_ptr<InputOLH>> olharray;
    std::vector<std::unique_ptr<InputHT>> lhtarray;
    std::vector<std::unique_ptr<InputRR>> uirrarray;
    std::stack<std::tuple < vector<int>, node_b*>> nodeStack;


    std::vector<std::vector<uint>>trueindexcache;



    std::vector <float> phi;


    std::default_random_engine generator;
    std::uniform_int_distribution<int> uni_level_sampler;
    //std::bernoulli_distribution bern_dist_rr;




};

#endif