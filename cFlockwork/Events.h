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
#ifndef __EVENTS__
#define __EVENTS__

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

using namespace std;

// ================================================ header =======================================

void rewire(
                 vector < set < size_t > * > & G, //Adjacency matrix
                 double Q,       //probability to connect with neighbors of neighbor
                 default_random_engine & generator, 
                 uniform_real_distribution<double> & distribution,
                 double & mean_degree,
                 set < pair < size_t, size_t > > & SI_E, //edge list of SI links
                 const vector < size_t > & node_status
           );

pair < vector < pair < size_t, size_t > >, vector < pair < size_t, size_t > > > 
    rewire_P(
                 vector < set < size_t > * > & G, //Adjacency matrix
                 double P,       //probability to connect with neighbors of neighbor
                 default_random_engine & generator, 
                 uniform_real_distribution<double> & distribution,
                 double & mean_degree,
                 set < pair < size_t, size_t > > & SI_E, //edge list of SI links
                 const vector < size_t > & node_status
           );

pair < vector < pair < size_t, size_t > >, vector < pair < size_t, size_t > > > 
    rewire_P_neighbor_affinity(
                 vector < set < size_t > * > & G, //Adjacency matrix
                 double P,       //probability to connect with neighbors of neighbor
                 vector < pair < vector < size_t >, vector < double > > > &neighbor_affinity,
                 vector < double > &total_affinity,
                 default_random_engine & generator, 
                 uniform_real_distribution<double> & distribution,
                 double & mean_degree,
                 set < pair < size_t, size_t > > & SI_E, //edge list of SI links
                 const vector < size_t > & node_status,
                 const bool use_preferential_node_selection = false
            );

void random_rewire(
                 vector < set < size_t > * > & G, //Adjacency matrix
                 default_random_engine & generator, 
                 uniform_real_distribution<double> & distribution,
                 set < pair < size_t, size_t > > & SI_E, //edge list of SI links
                 const vector < size_t > & node_status,
                 vector < size_t > & node_ints
           );

void infect(
                 vector < set < size_t > * > & G, //Adjacency matrix
                 default_random_engine & generator, 
                 uniform_real_distribution<double> & distribution,
                 set < pair < size_t, size_t > > & SI_E, //edge list of SI links
                 vector < size_t > & node_status,
                 set < size_t > & infected
           );

void SIS_recover(
                 vector < set < size_t > * > & G, //Adjacency matrix
                 default_random_engine & generator, 
                 uniform_real_distribution<double> & distribution,
                 set < pair < size_t, size_t > > & SI_E, //edge list of SI links
                 vector < size_t > & node_status,
                 set < size_t > & infected
           );

void SIR_recover(
                 vector < set < size_t > * > & G, //Adjacency matrix
                 default_random_engine & generator, 
                 uniform_real_distribution<double> & distribution,
                 set < pair < size_t, size_t > > & SI_E, //edge list of SI links
                 vector < size_t > & node_status,
                 set < size_t > & infected
           );

void SIRS_recover(
                 vector < set < size_t > * > & G, //Adjacency matrix
                 default_random_engine & generator, 
                 uniform_real_distribution<double> & distribution,
                 set < pair < size_t, size_t > > & SI_E, //edge list of SI links
                 vector < size_t > & node_status,
                 set < size_t > & infected,
                 set < size_t > & recovered
           );

void become_susceptible (
                 vector < set < size_t > * > & G, //Adjacency matrix
                 default_random_engine & generator, 
                 uniform_real_distribution<double> & distribution,
                 set < pair < size_t, size_t > > & SI_E, //edge list of SI links
                 vector < size_t > & node_status,
                 set < size_t > & recovered
           );

#endif
