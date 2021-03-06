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

#ifndef __MARKOV_DYNAMICS_H__
#define __MARKOV_DYNAMICS_H__

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
    markov_on_edge_lists(
            edge_lists & edg_lst,
            T & this_markov_object,
            double max_dt,
            bool is_static = false,
            bool verbose = false
            )
{
    if (verbose)
    {
        cout << "got markov object with properties" << endl;
        cout << "N = " << this_markov_object.N << endl;
        this_markov_object.print();
    }

    // reset the simulation object
    this_markov_object.reset();

    // initialize time variables
    double t0 = edg_lst.t[0];
    double tmax = edg_lst.tmax;
    double t_network_total = tmax - t0;
    double t = t0;
    double dt_to_next_change;

    this_markov_object.set_initial_time(t0);

    double t_simulation = this_markov_object.t_simulation;

    // initialize a graph
    size_t N = edg_lst.N;
    vector < set < size_t > > G(N);

    // create iterators
    auto next_edges = (edg_lst.edges).end();
    auto next_time = (edg_lst.t).end();
    long loop_count = -1;

    if (verbose)
    {
        cout << "initialized simulation with t-t0 = " << t-t0 << endl;
        cout << "empty graph of size = " << G.size() << endl;
        cout << "t_simulation = " << t_simulation << endl;
        cout << "simulation_ended = " << this_markov_object.simulation_ended() << endl;
    }

    while ( (t-t0 < t_simulation) and (not this_markov_object.simulation_ended()) )
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

            if ((not is_static) or (t == t0))
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
            if (not is_static)
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

        // update all the stuff that might have changed for this markov object
        // because G changed
        if ((not is_static) or (t == t0))
        {
            this_markov_object.update_network(G,t);
        }

        while ((dt_to_next_change > max_dt) and (t-t0+max_dt < t_simulation) and (not this_markov_object.simulation_ended()) )
        {
            this_markov_object.step(t, max_dt);
            dt_to_next_change -= max_dt;
            t += max_dt;
        }

        if ((t-t0+dt_to_next_change < t_simulation) and (not this_markov_object.simulation_ended()))
        {
            this_markov_object.step(t, dt_to_next_change);
        }
        t += dt_to_next_change;

        if (verbose)
        {
            cout << "the upcoming network time = " << this_next_time << endl;
            cout << "the upcoming network time with loops = " << this_next_time + t_network_total * loop_count << endl;
            cout << "dt_to_next_change = " << dt_to_next_change << endl;
            cout << "updated graph in markov instance" << endl;
            this_markov_object.print();
            cout << "========== PERFORMING markov =========" << endl;
        }

    }
}


template <typename T>
void 
    markov_on_edge_changes(
            edge_changes & edg_chg,
            T & this_markov_object,
            double max_dt,
            bool verbose = false
            )
{

    // reset the simulation object
    this_markov_object.reset();

    // initialize time variables
    double t0 = edg_chg.t0;
    double tmax = edg_chg.tmax;
    double t_network_total = tmax - t0;
    double t = t0;
    double dt_to_next_change;

    this_markov_object.set_initial_time(t0);

    double t_simulation = this_markov_object.t_simulation;

    // initialize a graph and pass it to the markov object
    size_t N = edg_chg.N;
    vector < set < size_t > > G(N);

    // create iterators
    auto next_edges_in = (edg_chg.edges_in).end();
    auto next_edges_out = (edg_chg.edges_out).end();
    auto next_time = (edg_chg.t).end();
    long loop_count = -1;

    while ( (t-t0 < t_simulation) and (not this_markov_object.simulation_ended()) )
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

            // update all the stuff that might have changed for this markov object
            // because G changed
            this_markov_object.update_network(G,t);
            this_markov_object.update_observables(t);
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

            // update all the stuff that might have changed for this markov object
            // because G changed
            this_markov_object.update_network(G,t);
            this_markov_object.update_observables(t);

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

        if (verbose)
        {
            cout << "the upcoming network time = " << this_next_time << endl;
            cout << "the upcoming network time with loops = " << this_next_time + t_network_total * loop_count << endl;
            cout << "dt_to_next_change = " << dt_to_next_change << endl;
            cout << "updated graph in markov instance" << endl;
            this_markov_object.print();
            cout << "========== PERFORMING markov =========" << endl;
        }

        while ((dt_to_next_change > max_dt) and (t-t0+max_dt < t_simulation) and (not this_markov_object.simulation_ended()) )
        {
            this_markov_object.step(t,max_dt);
            dt_to_next_change -= max_dt;
            t += max_dt;
        }

        if ((t-t0+dt_to_next_change < t_simulation) and (not this_markov_object.simulation_ended()))
        {
            this_markov_object.step(t, dt_to_next_change);
        }
        t += dt_to_next_change;
    }
}

#endif
