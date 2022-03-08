/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   node_b.h
 * Author: godzilla
 *
 * Created on 09 September 2018, 10:05
 */

#ifndef NODE_B_H
#define NODE_B_H

#include <vector>

class node_b {
public:

    float pert_freq_irr = 0.0f;


    float level_olh = 0.0f;
    float level_ht = 0.0f;

    uint hid=0;
    int id = 0;
    int b = 2;
    int level = 0;
    node_b *parent = NULL;
    std::vector <int> range = {0, 0};
    std::vector<node_b*> children;

    node_b(int l, int r, uint id, int b) {
        range[0] = l;
        range[1] = r;
        this->id = id;
        this->b = b;
        children = std::vector<node_b*>(b, NULL);
        //children.reserve(b);
    }

    node_b(uint b) {
        this->id = 0;
        this->b = b;
        children.reserve(b);
        children = std::vector<node_b*>(b, NULL);
    }

    void clearcounts() {

        pert_freq_irr = 0.0f;
        level_olh = 0.0f;
        level_ht = 0.0f;


    }
};


#endif /* NODE_B_H */

