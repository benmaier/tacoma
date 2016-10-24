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
#ifndef __EPIFLUTILS__
#define __EPIFLUTILS__

#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <set>
#include <utility>
#include <cmath>
#include <numeric>
#include <random>
#include <ctime>
#include <cstdlib>
#include <tuple>

using namespace std;

void choose (const size_t N, size_t &first, size_t &second, const double r1, const double r2);

pair <size_t,size_t> get_sorted_pair(size_t i, size_t j);

size_t arg_choose_from_vector(
        vector < double > const & weights, 
        default_random_engine & generator, 
        uniform_real_distribution<double> & distribution
      );

void get_gillespie_tau_and_event( 
                vector < double > const & rates,
                double & tau,
                size_t & event,
                default_random_engine & generator, 
                uniform_real_distribution<double> & distribution
             );

vector<size_t>::iterator choose_random_unique(
        vector<size_t>::iterator begin, 
        vector<size_t>::iterator end, 
        size_t num_random,
        default_random_engine & generator,
        uniform_real_distribution<double> & distribution
        ); 

void get_gillespie_tau_and_event_with_varying_gamma( 
                vector < double > & standard_rates,
                vector < pair < double, double > > const & gamma,
                double t0,
                double t_max,
                size_t & ti,
                double & tau,
                size_t & event,
                default_random_engine & generator, 
                uniform_real_distribution<double> & distribution
             );
#endif
