/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Timer.h
 * Author: godzilla
 *
 * Created on 08 September 2018, 21:04
 */
#pragma once
#ifndef TIMER_H
#define TIMER_H
#include <fstream>
#include <chrono>

class Timer {
    using clk = std::chrono::steady_clock;
    using microseconds = std::chrono::microseconds;

    clk::time_point tsb;
    clk::time_point tse;

public:

    void clear();
    void start();
    void stop();

    friend std::ostream& operator<<(std::ostream & o, const Timer & timer) {
        return o << timer.secs();
    }

    // return time difference in seconds

    double secs() const;
};






#endif /* TIMER_H */

