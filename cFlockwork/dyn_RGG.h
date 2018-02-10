/* 
 * The MIT License (MIT)
 * Copyright (c) 2018, Benjamin Maier
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

#ifndef __DYN_RGG_H__
#define __DYN_RGG_H__

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

void dyn_RGG_update_positions(vector < pair < double, double > > &pos,
                      default_random_engine &generator,
                      uniform_real_distribution<double> &uni_distribution,
                      const double &step_distance
                    );
        
void dyn_RGG_write_edge_list(
                   vector < pair < size_t, size_t > > & edges,
                   vector < pair < double, double > > const &pos,
                   const double &R,
                   const bool &PBC_distance
        );

size_t get_edge_int(pair < size_t, size_t > &edge, const size_t N);

edge_lists_with_histograms
     dynamic_RGG(
             const size_t N,       //number of nodes
             const size_t t_run_total,
             double step_distance,
             const double mean_link_duration,
             const double critical_density,
             const bool   PBC_distance,
             const bool   record_sizes_and_durations,
             const size_t seed,
             const bool verbose
        );

#endif
