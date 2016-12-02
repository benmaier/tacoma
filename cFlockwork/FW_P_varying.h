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

#ifndef __FW_P_VARYING_H__
#define __FW_P_VARYING_H__

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

edge_changes
     flockwork_P_varying_rates(
                 vector < pair < size_t, size_t > > E, //edgelist
                 const size_t N,       //number of nodes
                 vector < double > P,       //probability to reconnect after cutting
                 const double t_run_total,
                 vector < pair < double, double > > rewiring_rate,
                 const double tmax,
                 const bool   use_random_rewiring,
                 const bool   equilibrate_flockwork,
                 const size_t seed
        );

edge_changes
     flockwork_P_varying_rates_neighbor_affinity(
                 vector < pair < size_t, size_t > > E, //edgelist
                 const size_t N,       //number of nodes
                 vector < double > P,       //probability to reconnect after cutting
                 const double t_run_total,
                 vector < pair < double, double > > rewiring_rate,
                 vector < pair < vector < size_t >, vector < double > > > neighbor_affinity,
                 const double tmax,
                 const bool   use_random_rewiring,
                 const bool   equilibrate_flockwork,
                 const bool   use_preferential_node_selection,
                 const bool   use_unweighted_k_for_selection,
                 const size_t seed
        );


#endif
