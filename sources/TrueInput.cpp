
//#pragma once
#include <iostream>
#include "../Headers/helpers.h"
#include "../Headers/TrueInput.h"

void TrueInput::normalize()
{
        //norm_true_input=vector<float>(true_input.begin(), true_input.end());
        divBySum<std::vector<double>>(true_input);
        //cout << norm_true_input << endl;
}
void TrueInput::reset_object()
{
        for (size_t i = 0; i < D; i++)
        {
                true_input[i] = 0.0;
        }
        std::cout << "input object reset!" << std::endl;
}

