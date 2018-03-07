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

#include "SI.h"

using namespace std;
//namespace sis = SI;

void SI::update_network(
                    vector < set < size_t > > &_G,
                    double t
                    )
{
    // map this network to the current network
    G = &_G;

    // compute degree and R0
    size_t number_of_edges_times_two = 0;

    for(auto const &neighbors: *G)
        number_of_edges_times_two += neighbors.size();

    mean_degree = number_of_edges_times_two / (double) N;

    // update infected
    SI_edges.clear();

    // for each infected, check its neighbors.
    // if the neighbor is susceptible, push it to the
    // endangered susceptible vector
    for(auto const &inf: infected)
        for(auto const &neighbor_of_infected: (*G)[inf])
            if (node_status[neighbor_of_infected] == EPI::S)
                SI_edges.push_back( make_pair( inf, neighbor_of_infected ) );

    // update the arrays containing the observables
    update_observables(t);
}

void SI::get_rates_and_Lambda(
                    vector < double > &_rates,
                    double &_Lambda
                  )
{
    // delete current rates
    rates.clear();

    // compute rates of infection
    rates.push_back(infection_rate * SI_edges.size());

    // return those new rates
    _rates = rates;
    _Lambda = accumulate(rates.begin(),rates.end(),0.0);
}

void SI::make_event(
                size_t const &event,
                double t
               )
{
    if (event == 0)
        infection_event();
    else
        throw length_error("SI: chose event larger than rate vector which should not happen.");

    update_observables(t);
}

void SI::infection_event()
{
    // initialize uniform integer random distribution
    uniform_int_distribution < size_t > random_susceptible(0,SI_edges.size()-1);

    // find the index of the susceptible which will become infected
    size_t this_susceptible_index = random_susceptible(generator);

    // get the node number of this susceptible
    size_t this_susceptible = (SI_edges.begin() + this_susceptible_index)->second;

    // save this node as an infected
    infected.push_back(this_susceptible);

    // change node status of this node
    node_status[this_susceptible] = EPI::I;

    // erase all edges in the SI set where this susceptible is part of
    SI_edges.erase( 
            remove_if( 
                    SI_edges.begin(), 
                    SI_edges.end(),
                [&this_susceptible](const pair < size_t, size_t > & edge) { 
                        return edge.second == this_susceptible;
                }),
            SI_edges.end() 
        );

    // push the new SI edges
    size_t & this_infected = this_susceptible;
    for(auto const &neighbor: (*G)[this_infected])
        if (node_status[neighbor] == EPI::S)
            SI_edges.push_back( make_pair(this_infected, neighbor) );
}

void SI::update_observables(
                double t
               )
{
    // compute SI
    _SI.push_back(SI_edges.size());

    // compute I
    I.push_back(infected.size());

    // push back time
    time.push_back(t);
}

