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
#include "measurements.h"

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
#include <tuple>

using namespace std;

group_sizes_and_durations
     measure_group_sizes_and_durations(
             edge_lists &list_of_edge_lists,
             const bool verbose
        )
{
    // get references to edge_list and time
    vector < vector < pair < size_t, size_t > > > & all_edges = list_of_edge_lists.edges;
    vector < double > & time = list_of_edge_lists.t;
    size_t N = list_of_edge_lists.N;

    // create two graphs
    vector < set < size_t > > G(N);

    // create iteratos
    auto it_edge_lists = all_edges.begin();
    auto it_time = time.begin();
    double t0 = time[0];
    double tmax = list_of_edge_lists.tmax;
    double old_t;
    if (tmax < time.back())
        throw domain_error("The value tmax is smaller than the last value in the time list.");


    // container to track activities of nodes
    vector < double > last_time_active(N,t0);

    // maps and sets for measuring edge and group durations
    set < size_t > initial_edges;
    map < size_t, double > current_edges;
    auto current_edge_iterator = current_edges.begin();
    vector < double > contact_durations;
    vector < map < size_t, size_t > > size_histograms;
    vector < vector < double > > group_durations(N+1);

    // vectors and sets dealing with the change of components
    vector < set < size_t > > old_components;
    vector < set < size_t > > components;
    map < pair < size_t, size_t >, edge_weight > social_network;

    // iterate through all states
    while (it_edge_lists != all_edges.end() and it_time != time.end())
    {
        // get time of current state
        double t = *it_time;

        // get a sorted edge list
        vector < pair < size_t, size_t > > & these_edges = (*it_edge_lists);
        for(auto &edge: these_edges)
        {
            if (edge.first > edge.second)
                swap(edge.first,edge.second);
        }

        if (verbose) 
        {
            cout << " ============ " << endl;
            cout << "t = " << t << endl;
            cout << " creating Graph from " << these_edges.size() << " edges" << endl;
        }

        // get the current graph
        graph_from_edgelist(G,these_edges);

        if (verbose)
        {
            cout << " created Graph with degree sequence ";
            for( auto const &neighbors: G)
            {
                cout << neighbors.size() << " ";
            }
            cout << endl;
        }

        // get the current components and the size histogram
        map < size_t, size_t > this_size_histogram;
        components.clear();
        get_components_and_size_histogram(components,this_size_histogram,G);
        size_histograms.push_back(this_size_histogram);

        if (verbose)
        {
            cout << "found " << components.size() << " components" << endl;
        }

        if (t == t0) 
        {
            // initialize the measurement containers
            for(auto const &edge: these_edges)
            {
                size_t edge_int = get_edge_int(edge,N);
                initial_edges.insert(edge_int);
                current_edges[edge_int] = t;
            }
        }
        else
        {
            // =============== CHECK GROUP DIFFERENCES ================
            //
            
            //
            // for each node, get a pointer to its new component
            vector < set < size_t > * > new_component_of_node(N);
            for(auto &component: components)            
                for(auto const &node: component)
                    new_component_of_node[node] = &component;

            // for each component
            for(auto const &old_component: old_components)
            {
                // get the first node
                size_t first_node = *(old_component.begin());

                if (verbose)
                {
                    cout << "first node of component is node " << first_node << endl;
                }

                // get this node's new component
                set < size_t > * new_group = new_component_of_node[first_node];

                // if those components are not equal, the component changed
                if ( 
                        (old_component.size() != new_group->size()) or
                        (old_component != (*new_group) )
                   )
                {
                    // get the time difference to the last time a node in this
                    // component was active and push it to the group duration vector
                    // for the old component size.
                    // BUT only do this if the component is not still from the initialization
                    // because this will bias the duration histogram
                    double const &last_time = last_time_active[first_node];
                    if (last_time != t0)
                        group_durations[old_component.size()].push_back(t - last_time);

                    // all nodes of the old component are now in a new component and
                    // hence became active
                    for(auto const &node: old_component)
                        last_time_active[node] = t;
                }
            }

            // ================= CHECK LINK DURATION ==================
            //
            // check if any of the new edges are actually new
            // by checking wether or not they are in the current_edges dictionary
            for(auto const &edge: these_edges)
            {
                size_t edge_int = get_edge_int(edge,N);
                const bool already_in = current_edges.find(edge_int) != current_edges.end();
                if (not already_in)
                {
                    current_edges[edge_int] = t;
                }
            }

            // check if any of the old edges are not anymore in the new edges
            vector < size_t > edges_to_delete; 
            for(current_edge_iterator = current_edges.begin();
                current_edge_iterator != current_edges.end(); 
                current_edge_iterator++)
            {
                size_t edge_int = current_edge_iterator->first;
                size_t i = edge_int / N;
                size_t j = edge_int % N;
                pair < size_t, size_t > edge = make_pair(i,j);

                const bool not_in = find(these_edges.begin(),
                                         these_edges.end(),
                                         edge
                                         ) == these_edges.end();
                if (not_in)
                {
                    const bool is_initial_edge =    initial_edges.find(edge_int) 
                                                 != initial_edges.end();
                    if (initial_edges.size()==0 ||
                        not is_initial_edge
                       )
                    {
                        size_t duration = t - current_edge_iterator->second;
                        contact_durations.push_back(duration);
                    }
                    else if (is_initial_edge) {
                        initial_edges.erase(edge_int);
                    }
                    edges_to_delete.push_back(edge_int);
                }
            }
            for( auto const &edge_int: edges_to_delete)
            {
                current_edges.erase(edge_int);
            }
        }
        // add this state to the social network
        // if `edge` is not in `social_network` a weighted edge will be
        // constructed with .value = 0.0;
        if (t > t0)
            for(auto const &edge: these_edges)
                social_network[edge].value += t - old_t;

        // advance iterators
        it_edge_lists++;
        it_time++;
        old_t = t;

        // copy this graph and the components for comparisons in next time slice
        old_components = components;
    }

    // add the remaining edges to the social network
    for ( current_edge_iterator = current_edges.begin();
          current_edge_iterator != current_edges.end();
          ++current_edge_iterator )
    {
        size_t i = current_edge_iterator->first / N;
        size_t j = current_edge_iterator->first % N;
        double last_time = current_edge_iterator->second;
        social_network[make_pair(i,j)].value += tmax - last_time;
    }

    map < pair < size_t, size_t >, double > aggregated_network;
    for( auto soc_it = social_network.begin();
         soc_it != social_network.end();
         soc_it++
       )
    {
        aggregated_network[soc_it->first] = (soc_it->second).value;
    }

    group_sizes_and_durations result;
    result.contact_durations = contact_durations;
    result.size_histograms = size_histograms;
    result.group_durations = group_durations;
    result.aggregated_network = aggregated_network;

    return result;
    
}


group_sizes_and_durations
     measure_group_sizes_and_durations_for_edge_changes(
             edge_changes &list_of_edge_changes,
             const bool verbose
        )
{
    // get references to edge_list and time
    vector < vector < pair < size_t, size_t > > > & all_edges_in = list_of_edge_changes.edges_in;
    vector < vector < pair < size_t, size_t > > > & all_edges_out = list_of_edge_changes.edges_out;
    vector < double > & time = list_of_edge_changes.t;

    size_t N = list_of_edge_changes.N;
    // create graph
    vector < set < size_t > > G(N);
    graph_from_edgelist(G,list_of_edge_changes.edges_initial);
    double t0 = list_of_edge_changes.t0;
    double tmax = list_of_edge_changes.tmax;
    if (tmax < time.back())
        throw domain_error("The value tmax is smaller than the last time in the time list.");

    // create iterators
    auto it_edges_in = all_edges_in.begin();
    auto it_edges_out = all_edges_out.begin();
    auto it_time = time.begin();

    // container to track activities of nodes
    vector < double > last_time_active(N,t0);

    // maps and sets for measuring edge and group durations
    set < size_t > initial_edges;
    map < size_t, double > current_edges;
    auto current_edge_iterator = current_edges.begin();
    vector < double > contact_durations;
    vector < vector < double > > group_durations(N+1);

    // vectors and sets dealing with the change of components
    vector < set < size_t > > old_components;
    vector < set < size_t > > components;

    // initialize edges
    for(auto const &edge: list_of_edge_changes.edges_initial)
    {
        size_t edge_int = get_edge_int(edge,N);
        initial_edges.insert(edge_int);
        current_edges[edge_int] = t0;
    }

    // size_difference histograms
    vector < map < size_t, long > > size_histogram_differences;
    vector < map < size_t, size_t > > size_histograms;
    map < size_t, size_t > old_size_histogram;
    get_components_and_size_histogram(components,old_size_histogram,G);
    size_histograms.push_back(old_size_histogram);
    old_components = components;

    // for calculation of aggregated network
    map < pair < size_t, size_t >, edge_weight > social_network;

    // iterate through all changes
    while (it_edges_in != all_edges_in.end() and it_time != time.end() and it_edges_out != all_edges_out.end() )
    {
        // get time of current state
        double t = *it_time;

        // get a sorted edge list
        vector < pair < size_t, size_t > > & these_edges_in = (*it_edges_in);
        vector < pair < size_t, size_t > > & these_edges_out = (*it_edges_out);
        for(auto &edge: these_edges_out)
        {
            if (edge.first > edge.second)
                swap(edge.first,edge.second);
            size_t i = edge.first;
            size_t j = edge.second;
            G[i].erase(j);
            G[j].erase(i);
        }

        for(auto &edge: these_edges_in)
        {
            if (edge.first > edge.second)
                swap(edge.first,edge.second);
            size_t i = edge.first;
            size_t j = edge.second;
            G[i].insert(j);
            G[j].insert(i);
        }

        /*
        if (verbose) 
        {
            cout << " ============ " << endl;
            cout << "t = " << t << endl;
            cout << " creating Graph from " << these_edges.size() << " edges" << endl;
        }

        if (verbose)
        {
            cout << " created Graph with degree sequence ";
            for( auto const &neighbors: G)
            {
                cout << neighbors.size() << " ";
            }
            cout << endl;
        }
        */

        // get the current components and the size histogram
        map < size_t, size_t > this_size_histogram;
        components.clear();
        get_components_and_size_histogram(components,this_size_histogram,G);
        
        //get histogram difference because edge changes tend to be small so the histogram change should be small, too
        map < size_t, long > this_histogram_difference;
        for( auto const &this_hist_entry: this_size_histogram )
        {
            // get the size of the new group size
            size_t const & this_size = this_hist_entry.first;

            // if this entry is not in old_size_histogram, this means that there was a zero
            // and the difference is equal to the new count
            long difference = (long) this_hist_entry.second;

            // find this size in the old histogram
            auto old_size = old_size_histogram.find(this_size);
            if (old_size != old_size_histogram.end())
            {
                difference -= old_size->second;
            }
            this_histogram_difference[this_size] = difference;
        }

        //for the old histogram we only need to focus on entries which are not
        //in the new histogram since those were covered in the loop before
        for( auto const &old_hist_entry: old_size_histogram )
        {
            // get the size of the old group size
            size_t const & old_size = old_hist_entry.first;

            // find this size in the new histogram
            auto new_size = this_size_histogram.find(old_size);
            if (new_size == this_size_histogram.end())
            {
                // if entry does not exist, all groups of this
                // were deleted, so the difference is equal to
                // the old entry count
                long difference = -1 * ((long) old_hist_entry.second);
                this_histogram_difference[old_size] = difference;
            }
        }
        size_histogram_differences.push_back(this_histogram_difference);


        if (verbose)
        {
            cout << "found " << components.size() << " components" << endl;
        }

        // =============== CHECK GROUP DIFFERENCES ================
        //
        
        //
        // for each node, get a pointer to its new component
        vector < set < size_t > * > new_component_of_node(N);
        for(auto &component: components)            
            for(auto const &node: component)
                new_component_of_node[node] = &component;

        // for each component
        for(auto const &old_component: old_components)
        {
            // get the first node
            size_t first_node = *(old_component.begin());

            if (verbose)
            {
                cout << "first node of component is node " << first_node << endl;
            }

            // get this node's new component
            set < size_t > * new_group = new_component_of_node[first_node];

            // if those components are not equal, the component changed
            if ( 
                    (old_component.size() != new_group->size()) or
                    (old_component != (*new_group) )
               )
            {
                // get the time difference to the last time a node in this
                // component was active and push it to the group duration vector
                // for the old component size.
                // BUT only do this if the component is not still from the initialization
                // because this will bias the duration histogram
                double const &last_time = last_time_active[first_node];
                if (last_time != t0)
                    group_durations[old_component.size()].push_back(t - last_time);

                // all nodes of the old component are now in a new component and
                // hence became active
                for(auto const &node: old_component)
                    last_time_active[node] = t;
            }
        }

        // ================= CHECK LINK DURATION ==================
        //
        // check if any of the new edges are actually new
        // by checking wether or not they are in the current_edges dictionary
        // compute durations
        vector < size_t > edges_to_delete; 
        for(auto const &edge: these_edges_in)
        {
            size_t edge_int = get_edge_int(edge,N);
            current_edges[edge_int] = t;
            if (verbose)
                cout << "new_edge " << edge_int << " at time " << t << endl;
        }
        for(auto const &edge: these_edges_out)
        {
            size_t edge_int = get_edge_int(edge,N);
            const bool is_initial_edge =    initial_edges.find(edge_int) 
                                         != initial_edges.end();
            double duration = t - current_edges[edge_int];
            if (initial_edges.size() == 0 ||
                not is_initial_edge
               )
            {
                contact_durations.push_back(duration);
            }
            else if (is_initial_edge) {
                initial_edges.erase(edge_int);
            }

            // add this edge to the social network
            social_network[edge].value += duration;

            edges_to_delete.push_back(edge_int);
        }
        for( auto const &edge_int: edges_to_delete)
        {
            current_edges.erase(edge_int);
        }

        // advance iterators
        it_edges_in++;
        it_edges_out++;
        it_time++;

        // copy this graph and the components for comparisons in next time slice
        old_components = components;
        old_size_histogram = this_size_histogram;
    }

    // add the remaining edges to the social network
    for ( current_edge_iterator = current_edges.begin();
          current_edge_iterator != current_edges.end();
          ++current_edge_iterator )
    {
        size_t i = current_edge_iterator->first / N;
        size_t j = current_edge_iterator->first % N;
        double last_time = current_edge_iterator->second;
        social_network[make_pair(i,j)].value += tmax - last_time;
    }

    map < pair < size_t, size_t >, double > aggregated_network;
    for( auto soc_it = social_network.begin();
         soc_it != social_network.end();
         soc_it++
       )
    {
        aggregated_network[soc_it->first] = (soc_it->second).value;
    }


    group_sizes_and_durations result;
    result.contact_durations = contact_durations;
    result.size_histograms = size_histograms;
    result.size_histogram_differences = size_histogram_differences;
    result.group_durations = group_durations;
    result.aggregated_network = aggregated_network;

    return result;
    
}

