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

set < size_t > 
    get_edge_integer_set(
            vector < pair < size_t, size_t > > const & edges,
            size_t const & N
            )
{
    set < size_t > result;
    for(auto const &e: edges)
    {
        result.insert( hash_edge(e, N) );
    }

    return result;
}

edge_lists
     convert_edge_changes(
            edge_changes &list_of_edge_changes,
            const bool verbose
        )
{
    // get references to edge_list and time
    vector < vector < pair < size_t, size_t > > > & all_edges_in = list_of_edge_changes.edges_in;
    vector < vector < pair < size_t, size_t > > > & all_edges_out = list_of_edge_changes.edges_out;
    vector < double > & time = list_of_edge_changes.t;

    // create graph
    size_t N = list_of_edge_changes.N;
    vector < set < size_t > > G(N);
    graph_from_edgelist(G,list_of_edge_changes.edges_initial);

    // set initial and final time
    double t0 = list_of_edge_changes.t0;
    double tmax = list_of_edge_changes.tmax;

    // create iterators
    auto it_edges_in = all_edges_in.begin();
    auto it_edges_out = all_edges_out.begin();
    auto it_time = time.begin();

    vector < vector < pair < size_t, size_t > > > new_edges;
    vector < double > new_time;
    new_time.push_back(t0);
    
    vector < pair < size_t, size_t > > these_edges = list_of_edge_changes.edges_initial;
    new_edges.push_back(these_edges);

    // iterate through all changes
    while (it_edges_in != all_edges_in.end() and it_time != time.end() and it_edges_out != all_edges_out.end() )
    {
        // get time of current state
        double t = *it_time;

        new_time.push_back(t);

        // get a sorted edge list
        vector < pair < size_t, size_t > > & these_edges_in = (*it_edges_in);
        vector < pair < size_t, size_t > > & these_edges_out = (*it_edges_out);

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
        vector < pair < size_t, size_t > > edge_list;
        edgelist_from_graph(edge_list,G);

        new_edges.push_back(edge_list);

        // advance iterators
        it_edges_in++;
        it_edges_out++;
        it_time++;
    }

    edge_lists new_edge_lists;

    new_edge_lists.N = N;
    new_edge_lists.tmax = tmax;
    new_edge_lists.t = new_time;
    new_edge_lists.edges = new_edges;

    return new_edge_lists;
    
}

edge_changes
     convert_edge_lists(
            edge_lists &list_of_edge_lists,
            const bool verbose
        )
{
    // get references to edge_list and time
    vector < vector < pair < size_t, size_t > > > & all_edges = list_of_edge_lists.edges;
    vector < double > & time = list_of_edge_lists.t;

    // set initial and final time
    double t0 = time.front();
    double tmax = list_of_edge_lists.tmax;

    size_t N = list_of_edge_lists.N;

    // create iterators
    auto it_edges = all_edges.begin();
    auto it_time = time.begin();

    vector < double > new_time;
    vector < pair < size_t, size_t > > initial_edges = *it_edges;
    vector < vector < pair < size_t, size_t > > > edges_in;
    vector < vector < pair < size_t, size_t > > > edges_out;
    set < size_t > last_edge_integers = get_edge_integer_set(initial_edges, N);

    it_edges++;
    it_time++;

    // iterate through all changes
    while (it_edges != all_edges.end() and it_time != time.end() )
    {
        // get time of current state
        double t = *it_time;
        new_time.push_back(t);

        // get a sorted edge list
        set < size_t > edge_integers = get_edge_integer_set(*it_edges, N);

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

        edges_in.push_back( these_edges_in );
        edges_out.push_back( these_edges_out );

        it_edges++;
        it_time++;
        last_edge_integers = edge_integers;
    }

    edge_changes new_edge_changes;

    new_edge_changes.t = new_time;
    new_edge_changes.t0 = t0;
    new_edge_changes.tmax = tmax;
    new_edge_changes.N = N;
    new_edge_changes.edges_in = edges_in;
    new_edge_changes.edges_out = edges_out;
    new_edge_changes.edges_initial = initial_edges;
    
    return new_edge_changes;
}
