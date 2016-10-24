/* 
 * The MIT License (MIT)
 * Copyright (c) 2016, Benjamin Maier
 *
 * Permission is hereby granted, free of charge, to any person 
 * obtaining a copy of this software and associated documentation 
 * files (the "Software"), to deal in the Software without 
 * restriction, including without limitation the rights to use, 
 * copy, modify, merge, publish, distribute, sublicense, and/or 
 * sell copies of the Software, and to permit persons to whom the 
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall 
 * be included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NON-
 * INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS 
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN 
 * AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF 
 * OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS 
 * IN THE SOFTWARE.
 */

#ifndef __FLSIR_VARYING_H__
#define __FLSIR_VARYING_H__

#include "Events.h"
#include "Utilities.h"
#include "ResultClasses.h"
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <set>
#include <utility>
#include <random>
#include <cmath>
#include <numeric>
#include <random>
#include <ctime>
#include <tuple>

using namespace std;

SIR_result
     SIR_varying_rates(
                 vector < pair < size_t, size_t > > E, //edgelist
                 const size_t N,       //number of nodes
                 vector < double > Q,       //probability to connect with neighbors of neighbor
                 const double infection_rate,
                 const double recovery_rate,
                 vector < pair < double, double > > rewiring_rate,
                 const double tmax,
                 const size_t t_run_total,
                 const size_t number_of_vaccinated,
                 const size_t number_of_infected,
                 const bool   use_random_rewiring,
                 const bool   equilibrate_flockwork,
                 const size_t seed
        );

#endif
