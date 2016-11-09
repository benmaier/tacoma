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

#ifndef __EQFLOCKWORK_H__
#define __EQFLOCKWORK_H__

#include <iostream>
#include <vector>
#include <set>
#include <random>
#include <cmath>
#include <numeric>
#include <random>
#include <ctime>
#include <tuple>
#include "Utilities.h"

using namespace std;

// ================================================ header =======================================

vector < pair < size_t, size_t > > 
     equilibrate_edgelist_seed(
                 vector < pair < size_t, size_t >  > E, //edgelist
                 const size_t N,       //number of nodes
                 const double Q,       //probability to connect with neighbors of neighbor
                 const size_t seed = 0,
                 size_t t_max = 0,
                 const bool use_Q_as_P = false
                );

vector < pair < size_t, size_t > > 
     equilibrate_edgelist_generator(
                 vector < pair < size_t, size_t > > E, //edgelist
                 const size_t N,       //number of nodes
                 const double Q,       //probability to connect with neighbors of neighbor
                 default_random_engine & generator, 
                 uniform_real_distribution<double> & distribution,
                 size_t t_max,
                 const bool use_Q_as_P = false
                 );

void equilibrate_neighborset(
                 vector < set < size_t > * > &G, //edgelist
                 const double Q,       //probability to connect with neighbors of neighbor
                 default_random_engine & generator, 
                 uniform_real_distribution<double> & distribution,
                 size_t t_max = 0,
                 const bool use_Q_as_P = false
                );

void flockwork_timestep(
                 vector < set < size_t >  * > &G, //Adjacency matrix
                 double Q,       //probability to connect with neighbors of neighbor
                 default_random_engine & generator, 
                 uniform_real_distribution<double> & distribution
                );

void flockwork_P_timestep(
                 vector < set < size_t  > * > &G, //Adjacency matrix
                 double P,       //probability to connect with neighbors of neighbor
                 default_random_engine & generator, 
                 uniform_real_distribution<double> & distribution
            );

vector < pair < size_t, size_t > > 
     simulate_flockwork(
                 vector < pair < size_t, size_t > > E, //edgelist
                 const size_t N,       //number of nodes
                 const double Q,       //probability to connect with neighbors of neighbor
                 const size_t seed,
                 size_t num_timesteps
                 );

vector < pair < size_t, size_t > > 
     simulate_flockwork_P(
                 vector < pair < size_t, size_t > > E, //edgelist
                 const size_t N,       //number of nodes
                 const double P,       //probability to connect with neighbors of neighbor
                 const size_t seed,
                 size_t num_timesteps
                 );
#endif
