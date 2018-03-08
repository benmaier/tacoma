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
#include "conversion.h"
#include "social_trajectories.h"

using namespace std;

edge_changes
     concatenate_edge_changes(
            vector < edge_changes > & ec_lists,
            const bool verbose
        )
{
    edge_changes new_ec;

    auto ec_it = ec_lists.begin();

    // create graph
    size_t N = 0;

    // find maximum N
    for(; ec_it != ec_lists.end(); ++ec_it)
        if (ec_it->N > N)
            N = ec_it->N;

    // rewind
    ec_it = ec_lists.begin();
    double t0 = ec_it->t0;

    new_ec.N = N;
    new_ec.t0 = t0;
    new_ec.edges_initial = ec_it->edges_initial;


    vector < vector < pair < size_t, size_t > > > new_edges_in;
    vector < vector < pair < size_t, size_t > > > new_edges_out;

    vector < double > new_time;
    
    vector < pair < size_t, size_t > > these_edges = ec_it->edges_initial;

    double last_tmax = t0;

    // iterate through all edge_changes
    while ( ec_it != ec_lists.end() )
    {
        vector < set < size_t > > G(N);
        graph_from_edgelist(G,ec_it->edges_initial);

        // get time of current state

        // create iterators
        auto it_edges_in = (ec_it->edges_in).begin();
        auto it_edges_out = (ec_it->edges_out).begin();
        auto it_time = (ec_it->t).begin();

        while (     it_edges_in != (ec_it->edges_in).end() 
                and it_time != (ec_it->t).end() 
                and it_edges_out != (ec_it->edges_out).end() 
              )
        {

            // get a sorted edge list
            vector < pair < size_t, size_t > > these_edges_in = (*it_edges_in);
            vector < pair < size_t, size_t > > these_edges_out = (*it_edges_out);

            for(auto &edge: these_edges_out)
            {
                size_t i = edge.first;
                size_t j = edge.second;
                G[i].erase(j);
                G[j].erase(i);
            }

            for(auto &edge: these_edges_in)
            {
                size_t i = edge.first;
                size_t j = edge.second;
                G[i].insert(j);
                G[j].insert(i);
            }


            new_time.push_back( *it_time - (ec_it->t0) + last_tmax );
            new_edges_in.push_back(these_edges_in);
            new_edges_out.push_back(these_edges_out);
            
            // advance iterators
            it_edges_in++;
            it_edges_out++;
            it_time++;
        }

        // calculate edge changes between last state and next initial edge list
        if (ec_it + 1 != ec_lists.end())
        {
            vector < pair < size_t, size_t > > this_edge_list;
            edgelist_from_graph(this_edge_list,G);

            vector < pair < size_t, size_t > > next_edge_list = ec_it->edges_initial;

            // get a sorted edge list
            set < size_t > last_edge_integers = get_edge_integer_set(this_edge_list, N);
            set < size_t > edge_integers = get_edge_integer_set(next_edge_list, N);

            vector < size_t > incoming_edge_integers;
            vector < size_t > outgoing_edge_integers;

            // incoming edges are from the set (edges - last_edges) 
            set_difference(
                           edge_integers.begin(), edge_integers.end(),
                           last_edge_integers.begin(), last_edge_integers.end(),
                           back_inserter(incoming_edge_integers)
                          );

            // outgoing edges are from the set (last_edges - edges) 
            set_difference(
                           last_edge_integers.begin(), last_edge_integers.end(),
                           edge_integers.begin(), edge_integers.end(),
                           back_inserter(outgoing_edge_integers)
                          );

            vector < pair < size_t, size_t > > these_edges_in;
            vector < pair < size_t, size_t > > these_edges_out;

            for(auto const &e: incoming_edge_integers)
            {
                size_t j = e % N;
                size_t i = e / N;
                these_edges_in.push_back( make_pair(i, j) );
            }

            for(auto const &e: outgoing_edge_integers)
            {
                size_t j = e % N;
                size_t i = e / N;
                these_edges_out.push_back( make_pair(i, j) );
            }

            new_edges_in.push_back( these_edges_in );
            new_edges_out.push_back( these_edges_out );
            new_time.push_back(ec_it->tmax - ec_it->t0 + last_tmax );
        }

        last_tmax += ec_it->tmax - ec_it->t0;
        ec_it++;
    }

    double tmax = last_tmax;

    new_ec.tmax = tmax;
    new_ec.edges_in = new_edges_in;
    new_ec.edges_out = new_edges_out;
    new_ec.t = new_time;

    return new_ec;
}

edge_lists
     concatenate_edge_lists(
            vector < edge_lists > & el_lists,
            const bool verbose
            )
{
    edge_lists new_ec;

    auto ec_it = el_lists.begin();

    // create graph
    size_t N = 0;

    // find maximum N
    for(; ec_it != el_lists.end(); ++ec_it)
        if (ec_it->N > N)
            N = ec_it->N;

    // rewind
    ec_it = el_lists.begin();

    vector < vector < pair < size_t, size_t > > > new_edges;
    vector < double > new_time;
    
    double last_tmax = ec_it->t.front();

    // iterate through all edge_changes
    while ( ec_it != el_lists.end() )
    {
        // create iterators
        auto it_edges = (ec_it->edges).begin();
        auto it_time = (ec_it->t).begin();

        while (     
                    it_edges != (ec_it->edges).end() 
                and it_time != (ec_it->t).end() 
              )
        {
            vector < pair < size_t, size_t > > these_edges = *it_edges;

            new_edges.push_back(these_edges);
            new_time.push_back( *it_time - (ec_it->t.front()) + last_tmax );

            // advance iterators
            it_edges++;
            it_time++;
        }

        last_tmax += ec_it->tmax - ec_it->t.front();
        ec_it++;
    }

    double tmax = last_tmax;

    new_ec.N = N;
    new_ec.tmax = tmax;
    new_ec.edges = new_edges;
    new_ec.t = new_time;

    return new_ec;
}
