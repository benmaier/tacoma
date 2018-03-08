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

namespace EPI {
    enum {
        S = 0,
        I = 1,
        R = 2,
        E = 3
    };
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

struct edge_changes
{
    vector < double > t;
    vector < vector < pair < size_t, size_t > > > edges_out; 
    vector < vector < pair < size_t, size_t > > > edges_in;
    size_t N;
    vector < pair < size_t, size_t > > edges_initial; 
    double t0;
    double tmax;

    edge_changes(){};
    edge_changes(const edge_changes_with_histograms & other){
        copy_from(other);
    };

    void copy_from( const edge_changes_with_histograms & other) {

        edges_out = other.edges_out;
        edges_in = other.edges_in;
        edges_initial = other.edges_initial;

        t = other.t;
        N = other.N;
        t0 = other.t0;
        tmax = other.tmax;
    }
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

struct edge_lists
{
    vector < double > t;
    vector < vector < pair < size_t, size_t > > > edges;
    size_t N;
    double tmax;

    edge_lists(){};
    edge_lists(const edge_lists_with_histograms & other){
        copy_from(other);
    };

    void copy_from( const edge_lists_with_histograms & other) {

        edges = other.edges;

        t = other.t;
        N = other.N;
        tmax = other.tmax;
    }
};

struct group_sizes_and_durations
{

    vector < double > contact_durations;
    vector < map < size_t, size_t > > size_histograms;
    vector < map < size_t, long > > size_histogram_differences;
    vector < vector < double > > group_durations;
    vector < double > aggregated_size_histogram;
    map < pair < size_t, size_t >, double > aggregated_network;
};

struct edge_weight
{
    double value;
    edge_weight(): value(0.0){}
};

struct social_trajectory_entry {
    size_t hash;
    size_t size;
    vector < pair < double, double > > time_pairs;
};

struct edge_trajectory_entry {
    pair < size_t, size_t > edge;
    vector < pair < double, double > > time_pairs;
    double last_time_active;
    bool is_active;
};
#endif
