/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "../Headers/Timer.h"
void Timer::clear() {
    tsb = tse = clk::now();
}

void Timer::start() {
    tsb = clk::now();
}

void Timer::stop() {
    tse = clk::now();
}



// return time difference in seconds

double Timer::secs() const {
    if (tse <= tsb)
        return 0.0;
    auto d = std::chrono::duration_cast<microseconds>(tse - tsb);
    return d.count() / 1000000.0;
}
