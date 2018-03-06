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

#ifndef __DYN_GILLESPIE_H__
#define __DYN_GILLESPIE_H__

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

template <typename T>
void 
    gillespie_on_edge_lists(
            edge_lists & edg_lst,
            T & this_gillespie_object,
            bool verbose = false
            )
{
    if (verbose)
    {
        cout << "got gillespie object with properties" << endl;
        cout << "N = " << this_gillespie_object.N << endl;
        this_gillespie_object.print();
    }

    // deal with random numbers
    default_random_engine &generator = this_gillespie_object.generator;
    exponential_distribution<double> randexp(1.0);
    uniform_real_distribution<double> randuni(0.0,1.0);

    // initialize time variables
    double t0 = edg_lst.t[0];
    double tmax = edg_lst.tmax;
    double t_network_total = tmax - t0;
    double t = t0;
    double dt_to_next_change;

    double t_simulation = this_gillespie_object.t_simulation;

    // initialize a graph
    size_t N = edg_lst.N;
    vector < set < size_t > > G(N);

    // inititalize rate containers
    vector < double > rates;
    double Lambda;

    // create iterators
    auto next_edges = (edg_lst.edges).end();
    auto next_time = (edg_lst.t).end();
    long loop_count = -1;

    // draw time to first event
    double tau = randexp(generator);

    if (verbose)
    {
        cout << "initialized simulation with t-t0 = " << t-t0 << endl;
        cout << "empty graph of size = " << G.size() << endl;
        cout << "t_simulation = " << t_simulation << endl;
        cout << "simulation_ended = " << this_gillespie_object.simulation_ended() << endl;
    }

    while ( (t-t0 < t_simulation) and (not this_gillespie_object.simulation_ended()) )
    {

        if (verbose)
            cout << "=================== UPDATING THE GRAPH ================ " << endl;
        
        if (next_time == (edg_lst.t).end())
        {
            // calculate the current time
            // technically, this is unnecessary since we've covered
            // that in 
            //              t -= tau / Lambda
            // but this prevents the buildup of errors
            double this_time = tmax;
            t = this_time + t_network_total * loop_count;

            // if the upcoming edges is the end of the edges vector, this means
            // we've reached tmax. and at tmax, the network is getting looped
            // so we increase the loop count and generate G from the initial
            // edge list at edges.begin()
            loop_count += 1;
            graph_from_edgelist(G,*((edg_lst.edges).begin()));

            // furthermore, the upcoming change will be the second in the
            // vectors
            next_edges = (edg_lst.edges).begin() + 1;
            next_time = (edg_lst.t).begin() + 1;

            if (verbose)
            {
                cout << "initialized edge list from the beginning" << endl;
                cout << "graph has size " << G.size() << endl;
                cout << "the current time is t = " << t << endl;
            }
        }
        else
        {
            // calculate the current time
            // technically, this is unnecessary since we've covered
            // that in 
            //              t -= tau / Lambda
            // but this prevents the buildup of errors
            double this_time = *next_time;
            t = this_time + t_network_total * loop_count;

            // advance from the edges vector
            graph_from_edgelist(G,*next_edges);

            //advance to the upcoming change
            next_time++;
            next_edges++;
            
            if (verbose)
            {
                cout << "got graph from edgelist " << endl;
                cout << "the current time is t = " << t << endl;
            }
        }

        // if the upcoming time change points to the end of the edges vector,
        // the upcoming change will happen at tmax
        double this_next_time;
        if (next_time == (edg_lst.t).end())
            this_next_time = tmax;
        else
            this_next_time = *next_time;

        // compute the next dt while keeping in mind that the whole thing
        // might have looped already
        dt_to_next_change = this_next_time + t_network_total * loop_count - t;

        // update all the stuff that might have changed for this Gillespie object
        // because G changed
        this_gillespie_object.update_network(G,t);

        if (verbose)
        {
            cout << "the upcoming network time = " << this_next_time << endl;
            cout << "the upcoming network time with loops = " << this_next_time + t_network_total * loop_count << endl;
            cout << "dt_to_next_change = " << dt_to_next_change << endl;
            cout << "updated graph in Gillespie instance" << endl;
            this_gillespie_object.print();
            cout << "========== PERFORMING GILLESPIE =========" << endl;
        }

        // get the updated_rates
        this_gillespie_object.get_rates_and_Lambda(rates,Lambda);

        // ======================== PERFORM GILLESPIE ==========================
        
        if (tau >= Lambda * dt_to_next_change)
        {
            tau -= Lambda * dt_to_next_change;
            t += dt_to_next_change;
            
            if (verbose)
                cout << "no events happening in this time bin" << endl;
        }
        else
        {
            // the fraction of this time step that is left unitl
            // the next network change
            double xi = 1.0;

            // while it's not yet the time to do a network update
            while ( (tau < xi * Lambda * dt_to_next_change) and (t-t0 < t_simulation) and (not this_gillespie_object.simulation_ended()) ) 
            {
                // update the fraction that is left of this time step
                // till the next network change
                xi -= tau / (Lambda * dt_to_next_change);

                // advance time
                t +=  tau / Lambda;

                if (verbose)
                    cout << "advanced time to t = " << t << endl;

                if ((t - t0) < t_simulation)
                {
                    // choose event that takes place here
                    size_t event = 0;
                    auto this_rate = rates.begin();
                    double rProduct = randuni(generator) * Lambda;
                    double sum_event = 0.0;

                    while ( (this_rate != rates.end() ) and 
                            not ( (sum_event < rProduct) and (rProduct <= sum_event+(*this_rate)) ) )
                    {
                        sum_event += (*this_rate);
                        event++;
                        this_rate++;
                    }

                    if (verbose)
                        cout << "event happening in this time bin = " << event << endl;

                    // let the event take place
                    this_gillespie_object.make_event(event,t);

                    // get the updated rates and lambda
                    this_gillespie_object.get_rates_and_Lambda(rates,Lambda);

                    // get time till next event
                    tau = randexp(generator);
                }
            }

            // update tau because we'll advance to the next network
            // change
            tau -= xi * Lambda * dt_to_next_change;
            t += xi * dt_to_next_change;
        }
    }
}


template <typename T>
void 
    gillespie_on_edge_changes(
            edge_changes & edg_chg,
            T & this_gillespie_object,
            bool verbose = false
            )
{

    // deal with random numbers
    default_random_engine &generator = this_gillespie_object.generator;
    exponential_distribution<double> randexp(1.0);
    uniform_real_distribution<double> randuni(0.0,1.0);

    // initialize time variables
    double t0 = edg_chg.t0;
    double tmax = edg_chg.tmax;
    double t_network_total = tmax - t0;
    double t = t0;
    double dt_to_next_change;

    double t_simulation = this_gillespie_object.t_simulation;

    // initialize a graph and pass it to the gillespie object
    size_t N = edg_chg.N;
    vector < set < size_t > > G(N);

    // inititalize rate containers
    vector < double > rates;
    double Lambda;

    // create iterators
    auto next_edges_in = (edg_chg.edges_in).end();
    auto next_edges_out = (edg_chg.edges_out).end();
    auto next_time = (edg_chg.t).end();
    long loop_count = -1;

    // draw time to first event
    double tau = randexp(generator);

    while ( (t-t0 < t_simulation) and (not this_gillespie_object.simulation_ended()) )
    {

        // ======================== UPDATE THE GRAPH ==========================
        
        if (next_time == (edg_chg.t).end())
        {

            // calculate the current time
            // technically, this is unnecessary since we've covered
            // that in 
            //              t -= tau / Lambda
            // but this prevents the buildup of errors
            double this_time = tmax;
            t = this_time + t_network_total * loop_count;

            // if the upcoming change is the end of the change vector, this means
            // we've reached tmax. and at tmax, the network is getting looped
            // so we increase the loop count and generate G from the initial
            // edge list.
            loop_count += 1;
            graph_from_edgelist(G,edg_chg.edges_initial);

            // furthermore, the upcoming change will be the first in the
            // change vectors
            next_edges_in = (edg_chg.edges_in).begin();
            next_edges_out = (edg_chg.edges_out).begin();
            next_time = (edg_chg.t).begin();
        }
        else
        {
            // calculate the current time
            // technically, this is unnecessary since we've covered
            // that in 
            //              t -= tau / Lambda
            // but this prevents the buildup of errors
            double this_time = *next_time;
            t = this_time + t_network_total * loop_count;

            //update edges in
            for(auto & edge: *next_edges_in)
            {
                size_t &i = edge.first;
                size_t &j = edge.second;

                if (i==j)
                    throw domain_error("self loop detected.");

                if (i>j)
                    swap(i,j);

                G[i].insert(j);
                G[j].insert(i);
            }
            //update edges out
            for(auto & edge: *next_edges_out)
            {
                size_t &i = edge.first;
                size_t &j = edge.second;

                if (i==j)
                    throw domain_error("self loop detected.");

                if (i>j)
                    swap(i,j);

                G[i].erase(j);
                G[j].erase(i);
            }

            //advance to the upcoming change
            next_time++;
            next_edges_in++;
            next_edges_out++;
        }

        // if the upcoming time change points to the end of the change vector,
        // the upcoming change will happen at tmax
        double this_next_time;

        if (next_time == (edg_chg.t).end())
            this_next_time = tmax;
        else
            this_next_time = *next_time;

        // compute the next dt while keeping in mind that the whole thing
        // might have looped already
        dt_to_next_change = this_next_time + t_network_total * loop_count - t;

        // update all the stuff that might have changed for this Gillespie object
        // because G changed
        this_gillespie_object.update_network(G,t);

        // get the updated_rates
        this_gillespie_object.get_rates_and_Lambda(rates,Lambda);

        // ======================== PERFORM GILLESPIE ==========================
        
        if (tau >= Lambda * dt_to_next_change)
        {
            tau -= Lambda * dt_to_next_change;
            t += dt_to_next_change;
        }
        else
        {
            // the fraction of this time step that is left unitl
            // the next network change
            double xi = 1.0;

            // while it's not yet the time to do a network update
            while ( (tau < xi * Lambda * dt_to_next_change) and (t-t0 < t_simulation) and (not this_gillespie_object.simulation_ended()) ) 
            {
                // update the fraction that is left of this time step
                // till the next network change
                xi -= tau / (Lambda * dt_to_next_change);

                // advance time
                t +=  tau / Lambda;

                if ((t - t0) < t_simulation)
                {
                    // choose event that takes place here
                    size_t event = 0;
                    auto this_rate = rates.begin();
                    double rProduct = randuni(generator) * Lambda;
                    double sum_event = 0.0;
                    while ( (this_rate != rates.end() ) and 
                            not ( (sum_event < rProduct) and (rProduct <= sum_event+(*this_rate)) ) 
                          )
                    {
                        sum_event += (*this_rate);
                        event++;
                        this_rate++;
                    }

                    // let the event take place
                    this_gillespie_object.make_event(event,t);

                    // get the updated rates and lambda
                    this_gillespie_object.get_rates_and_Lambda(rates,Lambda);

                    // get time till next event
                    tau = randexp(generator);
                }
            }

            // update tau because we'll advance to the next network
            // change
            tau -= xi * Lambda * dt_to_next_change;
            t += xi * dt_to_next_change;
        }
    }
}

#endif
