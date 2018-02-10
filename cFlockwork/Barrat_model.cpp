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
        
edge_changes
     ZSBB_model(
             vector < pair < size_t, size_t > > E, //edgelist
             const size_t N,       //number of nodes
             const double lambda,
             const double b0,
             const double b1,
             const double t_run_simulation,
             const size_t seed,
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
    vector < set < size_t > * > G;

    //count number of edges in Graph
    size_t number_of_edges = 0;

    for(size_t node=0; node<N; node++)
    {
        G.push_back(new set < size_t >);
        last_time_active.push_back(t);
    }

    //loop through edge list and push neighbors
    for(auto edge: E)
    {
        //get nodes belonging to that edge
        size_t i = edge.first;
        size_t j = edge.second;

        //check if j is already a neighbor of i
        const bool already_counted = G[i]->find(j) != G[i]->end();

        //check if edge has been added already
        if (!(already_counted) && (i!=j)) 
        {
            G[ i ]->insert( j );
            G[ j ]->insert( i );
            number_of_edges++;
        }
    }

    // if a node has zero neighbors he/she is a loner
    for(size_t node=0; node<N; node++)
    {
        if (G[node]->size() == 0)
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
        const bool is_isolated = G[i]->size() == 0;

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

            const bool invite_loner = (uni_distribution(generator) < 1.0 - lambda);

            if ( ((is_isolated) || (invite_loner)) && (loners.size()>0) )
            {
                if (verbose)
                    cout << "will invite a loner" << endl;

                // choose a loner to invite for a group

                // calculate probabilities for the loners
                vector < double > probabilities;
                size_t loner_index;
                for(loner_index = 0; loner_index<loners.size(); ++loner_index)
                {
                    double tau_ = t - last_time_active[loners[loner_index]];
                    probabilities.push_back( p_n(b0,tau_,N) );
                }

                // choose a loner proportional to their activation probability
                do {
                    loner_index = arg_choose_from_vector(
                                        probabilities,
                                        generator,
                                        uni_distribution
                                        );
                } while(loners[loner_index] == i);

                // activate the loner 
                size_t j = loners[loner_index];
                last_time_active[j] = t;

                if (verbose)
                    cout << "invited loner " << j << endl;

                // add loner to all nodes of the group besides i
                for(auto neigh_i : *G[i] ){
                    last_time_active[neigh_i] = t;
                    G[neigh_i]->insert(j);
                    G[j]->insert(neigh_i);
                    in.push_back( get_sorted_pair(neigh_i, j));
                }

                // add loner to i
                G[i]->insert(j);
                G[j]->insert(i);
                in.push_back( get_sorted_pair(i,j) );

                // remove formerly isolated nodes from the loner pool
                if (is_isolated)
                    remove_2_from_vector(loners,i,j);
                else
                    remove_from_vector(loners,j);

            } 
            else // i is member of a group and does not invite anybody
            {
                if (verbose)
                    cout << "will leave group " << endl;

                // => leaves the group
                for(auto neigh_i : *G[i] ){

                    out.push_back( get_sorted_pair(i,neigh_i) );
                    if (verbose)
                        cout << "current_neighbor: " << neigh_i << endl;

                    last_time_active[neigh_i] = t;
                    G[neigh_i]->erase(i);

                    // if i was part of a pair, i's neighbor is now a loner, too
                    if (G[neigh_i]->size() == 0)
                        loners.push_back(neigh_i);
                }

                // delete all edges
                G[i]->clear();
                loners.push_back(i);

                if (verbose)
                    cout << "left group " << endl;

            }

            time.push_back(t);
            edges_in.push_back(in);
            edges_out.push_back(out);
        }


    }

    edge_changes result;

    result.t = time;
    result.edges_out = edges_out;
    result.edges_in = edges_in;

    return result;
}

