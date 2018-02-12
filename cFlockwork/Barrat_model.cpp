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

#include "Utilities.h"
#include "ResultClasses.h"
#include "Barrat_model.h"

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

double p_n(const double b, const double tau, const size_t N)
{
    return b / (1.0 + tau / N);
}
        
edge_changes_with_histograms
     ZSBB_model(
             vector < pair < size_t, size_t > > E, //edgelist
             const size_t N,       //number of nodes
             const double lambda,
             const double b0,
             const double b1,
             const double t_run_simulation,
             const size_t seed,
             const bool record_sizes_and_durations,
             const bool verbose
        )
{

    size_t t = 0;

    //initialize random generators
    default_random_engine generator(seed);
    uniform_real_distribution<double> uni_distribution(0.,1.);
    uniform_int_distribution<size_t> random_node(0,N-1);

    vector < size_t > last_time_active;
    vector < size_t > loners;

    //initialize Graph vector
    vector < set < size_t > > G(N);

    //count number of edges in Graph
    size_t number_of_edges = 0;

    for(size_t node=0; node<N; node++)
    {
        last_time_active.push_back(t);
    }

    //loop through edge list and push neighbors
    for(auto edge: E)
    {
        //get nodes belonging to that edge
        size_t i = edge.first;
        size_t j = edge.second;

        //check if j is already a neighbor of i
        const bool already_counted = G[i].find(j) != G[i].end();

        //check if edge has been added already
        if (!(already_counted) && (i!=j)) 
        {
            G[ i ].insert( j );
            G[ j ].insert( i );
            number_of_edges++;
        }
    }

    // get current size histogram
    vector < map < size_t, int > > group_changes;

    // add 
    vector < size_t > durations;
    set < size_t > initial_edges;
    map < size_t, size_t > current_edges;
    map < size_t, size_t >::iterator current_edge_iterator;
    for(auto const &edge: E)
    {
        size_t edge_int = get_edge_int(edge,N);
        initial_edges.insert(edge_int);
        current_edges[edge_int] = t;
    }

    // if a node has zero neighbors he/she is a loner
    for(size_t node=0; node<N; node++)
    {
        if (G[node].size() == 0)
        {
            loners.push_back(node);
        }
    }

    if (verbose)
        cout << "initiated" << endl;


    vector < vector < pair <size_t,size_t> > > edges_out;
    vector < vector < pair <size_t,size_t> > > edges_in;
    vector < double > time;



    for(t = 1; t<t_run_simulation; t++)
    {
        if (verbose)
        {
            cout << " ============== " << endl;
            cout << "t = " << t << endl;
        }

        // choose random node
        size_t i = random_node(generator);
        double b;
        double tau = t - last_time_active[i];
        const bool is_isolated = G[i].size() == 0;

        if (verbose)
        {
            cout << "chose node " << i << endl;
            cout << "tau = " << tau << endl;
            cout << "is_isolated = " << is_isolated << endl;
        }

        if (is_isolated) {
            b = b0;
        } else {
            b = b1;
        }

        // determine if node will become active
        if (uni_distribution(generator) < p_n(b,tau,N))
        {
            if (verbose)
            {
                cout << "is active " << endl;
            }

            last_time_active[i] = t;

            // if node will update
            vector < pair <size_t,size_t> > in;
            vector < pair <size_t,size_t> > out;
            map < size_t, int > this_group_change;

            const bool invite_loner = (uni_distribution(generator) < 1.0 - lambda);

            if ( ((is_isolated) || (invite_loner)) && (loners.size()>1) )
            {
                if (verbose)
                {
                    cout << "will invite a loner from loner vector of size " << loners.size() << endl;
                    cout << "      loner vector has content ";
                    for (auto const &loner: loners)
                    {
                        cout << loner << " ";
                    }
                    cout << endl;
                    cout << "      corresponding to probabilities " << endl;
                }

                // choose a loner to invite to i's group

                // calculate probabilities for the loners
                vector < double > probabilities;
                size_t i_loner_index;
                size_t loner_index = 0;
                for(auto const &loner: loners)
                {
                    // only consider loners which are not i
                    if (loner != i)
                    {
                        double tau_ = t - last_time_active[loner];
                        probabilities.push_back( p_n(b0,tau_,N) );
                    } else {
                        i_loner_index = loner_index;
                    }

                    loner_index++;
                }

                if (verbose)
                {
                    for(auto const &prob: probabilities)
                    {
                        cout << prob << " ";
                    }
                    cout << endl;
                }

                // choose a loner proportional to their activation probability
                loner_index = arg_choose_from_vector(
                                    probabilities,
                                    generator,
                                    uni_distribution
                                    );

                // in case node i was a loner, too we have to shift the chosen loner
                // in order to properly ignore i
                if ((is_isolated) && (loner_index >= i_loner_index))
                    loner_index += 1;

                if (verbose)
                {
                    cout << "chose loner node " << loners[loner_index] << " with loner index " << loner_index << " to invite to group " << endl;
                }

                // activate the loner 
                size_t j = loners[loner_index];
                last_time_active[j] = t;

                if (verbose)
                    cout << "invited loner " << j << endl;

                // add loner to all nodes of the group besides i
                for(auto const &neigh_i : G[i] ){
                    last_time_active[neigh_i] = t;
                    G[neigh_i].insert(j);
                    G[j].insert(neigh_i);
                    in.push_back( get_sorted_pair(neigh_i, j));
                }

                // add loner to i
                G[i].insert(j);
                G[j].insert(i);
                in.push_back( get_sorted_pair(i,j) );

                // remove formerly isolated nodes from the loner pool
                if (is_isolated)
                {
                    remove_2_from_vector(loners,i,j);
                    this_group_change[1] = -2;
                    this_group_change[2] = +1;
                }
                else
                {
                    remove_from_vector(loners,j);
                    this_group_change[1] = -1;
                    this_group_change[G[i].size()] = -1;
                    this_group_change[G[i].size()+1] = +1;
                }

            } 
            else if (not is_isolated)// i is member of a group and does not invite anybody
            {
                if (verbose)
                    cout << "will leave group " << endl;

                bool const is_pair = G[i].size() == 1;
                size_t const old_group_size = G[i].size() + 1;

                // => leaves the group
                for(auto const &neigh_i : G[i] ){

                    out.push_back( get_sorted_pair(i,neigh_i) );
                    if (verbose)
                        cout << "current_neighbor: " << neigh_i << endl;

                    last_time_active[neigh_i] = t;
                    G[neigh_i].erase(i);

                    // if i was part of a pair, i's neighbor is now a loner, too
                    if (G[neigh_i].size() == 0)
                        loners.push_back(neigh_i);
                }

                if (is_pair)
                {
                    this_group_change[2] = -1;
                    this_group_change[1] = +2;
                }
                else
                {
                    this_group_change[old_group_size] = -1;
                    this_group_change[old_group_size-1] = +1;
                    this_group_change[1] = +1;
                }

                // delete all edges
                G[i].clear();
                loners.push_back(i);

                if (verbose)
                    cout << "left group " << endl;

            }

            time.push_back(t);
            edges_in.push_back(in);
            edges_out.push_back(out);

            if(record_sizes_and_durations)
            {
                // compute durations
                vector < size_t > edges_to_delete; 
                for(auto const &edge: in)
                {
                    size_t edge_int = get_edge_int(edge,N);
                    current_edges[edge_int] = t;
                    if (verbose)
                        cout << "new_edge " << edge_int << " at time " << t << endl;
                }
                for(auto const &edge: out)
                {
                    size_t edge_int = get_edge_int(edge,N);
                    const bool is_initial_edge =    initial_edges.find(edge_int) 
                                                 != initial_edges.end();
                    if (initial_edges.size()==0 ||
                        not is_initial_edge
                       )
                    {
                        size_t duration = t - current_edges[edge_int];
                        durations.push_back(duration);
                    }
                    else if (is_initial_edge) {
                        initial_edges.erase(edge_int);
                    }
                    edges_to_delete.push_back(edge_int);
                }
                for( auto const &edge_int: edges_to_delete)
                {
                    current_edges.erase(edge_int);
                }

                //compute histogram
                group_changes.push_back(this_group_change);
            }
        }


    }

    edge_changes_with_histograms result;

    result.t = time;
    result.edges_out = edges_out;
    result.edges_in = edges_in;
    result.group_changes = group_changes;
    result.durations = durations;

    return result;
}

