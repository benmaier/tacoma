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
#include <thread>
#include <map>

using namespace std;

size_t get_edge_int(pair < size_t, size_t > const &edge, const size_t N);

void choose (const size_t N, size_t &first, size_t &second, const double r1, const double r2);

pair <size_t,size_t> get_sorted_pair(size_t i, size_t j);

size_t arg_choose_from_vector(
        vector < double > const & weights, 
        mt19937_64 & generator, 
        uniform_real_distribution<double> & distribution
      );

size_t arg_choose_from_vector(
        vector < size_t > const & weights, 
        mt19937_64 & generator, 
        uniform_real_distribution<double> & distribution
      );

void get_gillespie_tau_and_event( 
                vector < double > const & rates,
                double & tau,
                size_t & event,
                mt19937_64 & generator, 
                uniform_real_distribution<double> & distribution
             );

vector<size_t>::iterator choose_random_unique(
        vector<size_t>::iterator begin, 
        vector<size_t>::iterator end, 
        size_t num_random,
        mt19937_64 & generator,
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
                mt19937_64 & generator, 
                uniform_real_distribution<double> & distribution
             );

void get_gillespie_tau_and_event_with_varying_gamma_for_each_node( 
                vector < double > & standard_rates,
                vector < pair < double, double > > const & gamma,
                vector < pair < double, vector < double > > > const & gamma_single_nodes,
                double t0,
                double t_max,
                size_t & i_t,
                double & tau,
                size_t & event,
                mt19937_64 & generator, 
                uniform_real_distribution<double> & distribution
             );

void remove_from_vector(vector <size_t> &vec, const size_t to_be_removed);
void remove_2_from_vector(vector <size_t> &vec, const size_t first_to_be_removed, const size_t second_to_be_removed);
/*
void randomly_seed_engine(
        mt19937_64 &generator
        );
*/
set < size_t > get_component_of_node(
            size_t &node,
            const vector < set < size_t > > &G
        );
void get_components_and_size_histogram(
        vector < set <size_t> > &components,
        map < size_t, size_t > &counter,
        const vector < set < size_t > > &G
        );
void get_component_size_histogram_from_edgelist(
        size_t N,
        vector < pair < size_t,size_t > > const &edge_list,
        map < size_t,size_t > &counter
        );
void get_component_size_histogram(
        map < size_t, size_t > &counter,
        const vector < set < size_t > > &G
        );
void add_nodes_belonging_to_this_component(
        size_t start_node,
        const vector < set < size_t > > &G,
        set < size_t > &comp,
        vector < bool > &already_visited
       );

void graph_from_edgelist(vector < set < size_t > > &G,
                         vector < pair < size_t, size_t > > &edge_list
                         );

void graph_from_edgelist(vector < set < size_t > * > &G,
                         vector < pair < size_t, size_t > > &edge_list
                         );

void edgelist_from_graph(
                         vector < pair < size_t, size_t > > &edge_list,
                         vector < set < size_t > > &G
                         );


void randomly_seed_engine(
        mt19937_64 &generator
        );

void seed_engine(
        mt19937_64 &generator, 
        size_t seed
        );

vector < set < size_t > > get_random_graph(
        size_t n,
        double p,
        mt19937_64 &generator
        );
#endif
