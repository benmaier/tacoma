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

#ifndef __SOC_TRAJ_H__
#define __SOC_TRAJ_H__

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
#include <functional>

using namespace std;

size_t get_group_integer(
            size_t const & N,
            set < size_t > const &component,
            map < size_t, size_t > &hash_to_int,
            size_t &this_hash
        );

vector < set < size_t > >
     binned_social_trajectory_from_edge_lists(
             edge_lists &list_of_edge_lists,
             size_t node,
             double dt = 0.0,
             size_t N_time_steps = 0,
             const bool verbose = false
        );

vector < set < size_t > >
     binned_social_trajectory_from_edge_changes(
             edge_changes &list_of_edge_changes,
             size_t node,
             double dt = 0.0,
             size_t N_time_steps = 0,
             const bool verbose = false
        );

vector < social_trajectory_entry >
     social_trajectory_from_edge_changes(
            edge_changes &list_of_edge_changes,
            size_t node,
            const bool verbose
        );

vector < social_trajectory_entry >
     social_trajectory_from_edge_lists(
            edge_lists &list_of_edge_lists,
            size_t node,
            const bool verbose
        );

size_t hash_edge(const pair < size_t, size_t > &p, const size_t &N);

size_t get_edge_integer(
            const size_t &N,
            pair < size_t, size_t > const &edge,
            map < size_t, size_t > &hash_to_int
        );

vector < edge_trajectory_entry >
        edge_trajectories_from_edge_lists(
            edge_lists &list_of_edge_lists,
            const bool verbose
        );

vector < edge_trajectory_entry >
        edge_trajectories_from_edge_changes(
            edge_changes &list_of_edge_changes,
            const bool verbose
        );
#endif
