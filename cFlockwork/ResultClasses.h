/* 
 * The MIT License (MIT)
 * Copyright (c) 2016-2018, Benjamin Maier
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

#ifndef __RESULTCLASSES_H__
#define __RESULTCLASSES_H__

#include "Events.h"
#include "Utilities.h"
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

struct SIR_result 
{
    vector < pair < double, size_t > > I_of_t; 
    vector < pair < double, size_t > > R_of_t; 
    vector < pair < double, size_t > > SI_of_t; 
    vector < pair < double, double > > R0_of_t;
    set < pair < size_t, size_t > > edge_list;
};

struct SIS_result 
{
    vector < pair < double, size_t > > I_of_t; 
    vector < pair < double, size_t > > SI_of_t; 
    vector < pair < double, double > > R0_of_t;
    set < pair < size_t, size_t > > edge_list;
};

struct edge_changes
{
    vector < double > t;
    vector < vector < pair < size_t, size_t > > > edges_out; 
    vector < vector < pair < size_t, size_t > > > edges_in;
    size_t N;
    vector < pair < size_t, size_t > > edges_initial; 
    double t0;
    double tmax;
};

struct edge_changes_with_histograms
{
    vector < double > t;
    vector < vector < pair < size_t, size_t > > > edges_out; 
    vector < vector < pair < size_t, size_t > > > edges_in;
    map < size_t, size_t > initial_size_histogram;
    vector < map < size_t, int > > group_changes;
    map < size_t, size_t > final_size_histogram;
    vector < size_t > contact_durations;
    vector < size_t > inter_contact_durations;
    vector < vector < size_t > > group_durations;

    size_t N;
    vector < pair < size_t, size_t > > edges_initial; 
    double t0;
    double tmax;
};

struct edge_lists
{
    vector < double > t;
    vector < vector < pair < size_t, size_t > > > edges;
    size_t N;
    double tmax;
};

struct edge_lists_with_histograms
{
    vector < double > t;
    vector < vector < pair < size_t, size_t > > > edges;
    vector < map < size_t, size_t > > size_histograms;
    vector < size_t > group_durations;
    size_t N;
    double tmax;
};

struct group_sizes_and_durations
{

    vector < double > contact_durations;
    vector < map < size_t, size_t > > size_histograms;
    vector < vector < double > > group_durations;
    map < pair < size_t, size_t >, double > aggregated_network;
};

struct edge_weight
{
    double value;
    edge_weight(): value(0.0){}
};


#endif
