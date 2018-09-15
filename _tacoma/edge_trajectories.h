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

#ifndef __EDGE_TRAJ_H__
#define __EDGE_TRAJ_H__

#include "Utilities.h"
#include "ResultClasses.h"
#include "social_trajectories.h"

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

struct event_map_entry {
    vector < pair < size_t, size_t > > edges_in;
    vector < pair < size_t, size_t > > edges_out;
};

edge_trajectories
        edge_trajectories_from_edge_lists(
            edge_lists &list_of_edge_lists,
            const bool return_edge_similarities = false,
            const bool verbose = false
        );

edge_trajectories
        edge_trajectories_from_edge_changes(
            edge_changes &list_of_edge_changes,
            const bool return_edge_similarities = false,
            const bool verbose = false
        );

edge_changes
    edge_trajectories_to_edge_changes(
                edge_trajectories &traj
            );
#endif
