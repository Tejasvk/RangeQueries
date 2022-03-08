
/* 
 * File:   node_b_l.h
 * Author: godzilla
 *
 * Created on 09 September 2018, 09:52
 */
#pragma once
#ifndef NODE_B_L_H
#define NODE_B_L_H
#include <iostream>
#include <vector>

class node_b_l {
public:
    double pert_freq_irr = 0.0;
    double level_ht = 0.0;    
    double level_olh = 0.0;
     
    uint hid = 0;
    int id = 0;
    int b = 2;
    int level = 0;
    node_b_l *parent = NULL;
    std::vector <int> range = {0, 0};
    std::vector<node_b_l*> children;

    node_b_l(int l, int r, uint id, int b) {
        range = {l, r};
        this->id = id;
        this->b = b;
        children = std::vector<node_b_l*>(b, NULL);
        //children.reserve(b);
    }

    node_b_l(uint b) {
        this->id = 0;
        this->b = b;
        //children.reserve(b);
        children = std::vector<node_b_l*>(b, NULL);
    }

    void clearcounts() {

        pert_freq_irr = 0.0;
        level_ht = 0.0;
        level_olh = 0.0;      

    }
};

#endif /* NODE_B_L_H */

