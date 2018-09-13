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
#include "slice.h"

using namespace std;

edge_lists
     slice_edge_lists(
             edge_lists &list_of_edge_lists,
             const double new_t0,
             const double new_tmax,
             const bool verbose
        )
{

    // get references to edge_list and time
    vector < vector < pair < size_t, size_t > > > & all_edges = list_of_edge_lists.edges;
    vector < double > &time = list_of_edge_lists.t;
    double tmax = list_of_edge_lists.tmax;
    double t0 = time[0];


    if (new_tmax > tmax)
        throw domain_error("The value new_tmax is greater than edge_lists.tmax.");

    if (new_t0 < t0)
        throw domain_error("The value new_t0 is lower than edge_lists.t[0].");

    if (new_t0 >= new_tmax)
        throw domain_error("new_t0 does not fullfil new_t0 < new_tmax");

    size_t N = list_of_edge_lists.N;

    // create iteratos
    auto it_edge_lists = all_edges.begin();
    auto it_time = time.begin();

    vector < vector < pair < size_t, size_t > > > new_list_of_edge_lists;
    vector < double > new_time;

    while( (it_time+1) != time.end() and *(it_time+1) <= new_t0 )
    {
        ++it_time;
        ++it_edge_lists;
    }

    new_list_of_edge_lists.push_back( *it_edge_lists );
    new_time.push_back( new_t0 );

    ++it_time;
    ++it_edge_lists;


    while( it_time != time.end() and *it_time < new_tmax )
    {
        new_list_of_edge_lists.push_back( *it_edge_lists );
        new_time.push_back( *it_time );
        ++it_time;
        ++it_edge_lists;
    }


    edge_lists result;
    result.tmax = new_tmax;
    result.N = N;
    result.edges = new_list_of_edge_lists;
    result.t = new_time;
    result.time_unit = list_of_edge_lists.time_unit;
    result.notes = list_of_edge_lists.notes;
    result.int_to_node = list_of_edge_lists.int_to_node;

    return result;
}


edge_changes
     slice_edge_changes(
             edge_changes &list_of_edge_changes,
             const double new_t0,
             const double new_tmax,
             const bool verbose
             )
{
    // get references to edge_list and time
    vector < vector < pair < size_t, size_t > > > & all_edges_in = list_of_edge_changes.edges_in;
    vector < vector < pair < size_t, size_t > > > & all_edges_out = list_of_edge_changes.edges_out;
    vector < double > &time = list_of_edge_changes.t;

    // set initial and final time
    double t0 = list_of_edge_changes.t0;
    double tmax = list_of_edge_changes.tmax;

    // create graph
    size_t N = list_of_edge_changes.N;
    vector < set < size_t > > G(N);
    graph_from_edgelist(G,list_of_edge_changes.edges_initial);

    // check times
    if (new_tmax > tmax)
        throw domain_error("The value new_tmax is greater than edge_lists.tmax.");

    if (new_t0 < t0)
        throw domain_error("The value new_t0 is lower than edge_lists.t[0].");

    if (new_t0 >= new_tmax)
        throw domain_error("new_t0 does not fullfil new_t0 < new_tmax");

    // create iteratos
    auto it_edges_in = all_edges_in.begin();
    auto it_edges_out = all_edges_out.begin();
    auto it_time = time.begin();

    vector < vector < pair < size_t, size_t > > > new_list_of_edge_lists;
    vector < double > new_time;

    while (it_time != time.end() and *(it_time) <= new_t0)
    {
        for(auto const &edge: *it_edges_in)
        {

            size_t const & i = edge.first;
            size_t const & j = edge.second;
            G[ i ].insert( j );
            G[ j ].insert( i );
        }

        for(auto const &edge: *it_edges_out)
        {
            size_t const & i = edge.first;
            size_t const & j = edge.second;
            G[ i ].erase( j );
            G[ j ].erase( i );
        }
        
        ++it_time;
        ++it_edges_in;
        ++it_edges_out;
    }

    //if (*it_time == new_t0)
    //{
    //    ++it_time;
    //    ++it_edges_in;
    //    ++it_edges_out;
    //}


    vector < pair < size_t, size_t > > new_initial_edge_list;
    vector < vector < pair < size_t, size_t > > > new_edges_in;
    vector < vector < pair < size_t, size_t > > > new_edges_out;

    edgelist_from_graph(new_initial_edge_list,G);

    while( it_time != time.end() and *it_time < new_tmax )
    {
        new_edges_in.push_back(*it_edges_in);
        new_edges_out.push_back(*it_edges_out);
        new_time.push_back( *it_time );

        ++it_time;
        ++it_edges_in;
        ++it_edges_out;
    }


    edge_changes result;
    result.tmax = new_tmax;
    result.t0 = new_t0;
    result.N = N;
    result.edges_initial = new_initial_edge_list;
    result.edges_in = new_edges_in;
    result.edges_out = new_edges_out;
    result.t = new_time;
    result.time_unit = list_of_edge_changes.time_unit;
    result.notes = list_of_edge_changes.notes;
    result.int_to_node = list_of_edge_changes.int_to_node;

    return result;
}
