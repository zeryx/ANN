/* 
 * File:   main.cpp
 * Author: James
 *
 * Created on January 20, 2014, 12:56 PM
 */

#include <boost/multi_array.hpp>
#include <cassert>
#include <cmath>
#include <random>
#include "structs.h"
#include "classdef.h"
#include "members.h"

double Activation_sig(double X);

int main() {
    Neural test(4, 1, 1, 100, 0.0001, 1000, 10^-6);
    test.prepare();
    test.callgen();
    return 0;
}

