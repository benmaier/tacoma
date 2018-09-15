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
#include "edge_trajectories.h"
#include "conversion.h"

using namespace std;

edge_trajectories
        edge_trajectories_from_edge_changes(
            edge_changes &list_of_edge_changes,
            const bool return_edge_similarities,
            const bool verbose
        )
{
    if (return_edge_similarities)
    {
        edge_lists el = convert_edge_changes(list_of_edge_changes);
        return edge_trajectories_from_edge_lists(el);
    }

    // get references to edge_list and time
    vector < vector < pair < size_t, size_t > > > & all_edges_in = list_of_edge_changes.edges_in;
    vector < vector < pair < size_t, size_t > > > & all_edges_out = list_of_edge_changes.edges_out;
    vector < double > & time = list_of_edge_changes.t;

    // create graph
    size_t N = list_of_edge_changes.N;

    // set initial and final time
    double t0 = list_of_edge_changes.t0;
    double tmax = list_of_edge_changes.tmax;

    if (verbose)
    {
        cout << "last time in array = " << time.back() << endl;
        cout << "              tmax = " << tmax << endl;
    }

    if (tmax < time.back())
        throw domain_error("The value tmax is smaller than the last time in the time list.");

    // create iterators
    auto it_edges_in = all_edges_in.begin();
    auto it_edges_out = all_edges_out.begin();
    auto it_time = time.begin();


    // vectors and sets dealing with the change of edges
    map < size_t, size_t > hash_to_int;

    vector < edge_trajectory_entry > trajectories;

    for(auto &edge: list_of_edge_changes.edges_initial)
    {
        if (edge.first > edge.second)
            swap(edge.first,edge.second);
        size_t edge_int = get_edge_integer(N,edge,hash_to_int);
        edge_trajectory_entry this_entry;
        if (edge_int == trajectories.size())
        {
            this_entry.edge = edge;
            this_entry.last_time_active = t0;
        }
        else
            throw domain_error("Doubling of edges in edge_changes.edges_initial");

        trajectories.push_back(this_entry);
    }


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

            size_t edge_int = get_edge_integer(N,edge,hash_to_int);

            if (edge_int == trajectories.size())
                throw length_error("an edge was in edges_out, but did not exist before");

            edge_trajectory_entry &this_entry = trajectories[edge_int];
            this_entry.time_pairs.push_back(make_pair(this_entry.last_time_active, t));
            this_entry.last_time_active = t0 - 1000.0;
        }

        for(auto &edge: these_edges_in)
        {
            if (edge.first > edge.second)
                swap(edge.first,edge.second);

            size_t edge_int = get_edge_integer(N,edge,hash_to_int);

            if (edge_int == trajectories.size())
            {
                edge_trajectory_entry this_entry;
                this_entry.edge = edge;
                this_entry.last_time_active = t;
                trajectories.push_back(this_entry);
            }
            else
            {
                trajectories[edge_int].last_time_active = t;
            }

        }

        it_edges_in++;
        it_edges_out++;
        it_time++;
    }

    for(auto &this_entry: trajectories)
    {
        if (this_entry.last_time_active >= t0)
            this_entry.time_pairs.push_back(make_pair(this_entry.last_time_active,tmax));
    }

    edge_trajectories traj;

    vector < tuple < size_t, size_t, double > > sim;

    traj.trajectories = trajectories;
    traj.edge_similarities = sim;

    traj.N = list_of_edge_changes.N;
    traj.t0 = list_of_edge_changes.t0;
    traj.tmax = list_of_edge_changes.tmax;
    traj.time_unit = list_of_edge_changes.time_unit;
    traj.notes = list_of_edge_changes.notes;
    traj.int_to_node = list_of_edge_changes.int_to_node;

    return traj;
    
}

edge_trajectories
        edge_trajectories_from_edge_lists(
            edge_lists &list_of_edge_lists,
            const bool return_edge_similarities,
            const bool verbose
        )
{
    //cout << "return_edge_similarities = " << return_edge_similarities << endl;
    //cout << "verbose = " << verbose << endl;

    // get references to edge_list and time
    vector < vector < pair < size_t, size_t > > > & all_edges = list_of_edge_lists.edges;
    vector < double > & time = list_of_edge_lists.t;

    // create graph
    size_t N = list_of_edge_lists.N;

    // set initial and final time
    double tmax = list_of_edge_lists.tmax;

    if (verbose)
    {
        cout << "last time in array = " << time.back() << endl;
        cout << "              tmax = " << tmax << endl;
    }

    if (tmax < time.back())
        throw domain_error("The value tmax is smaller than the last time in the time list.");

    // create iterators
    auto it_edges = all_edges.begin();
    auto it_time = time.begin();

    // vectors and sets dealing with the change of edges
    map < size_t, size_t > hash_to_int;

    vector < edge_trajectory_entry > these_edge_trajectories;

    map < pair < size_t, size_t >, edge_weight > similarities;

    // iterate through all changes
    while (it_edges != all_edges.end() and it_time != time.end() )
    {
        // get time of current and next state
        double t = *it_time;
        double next_time;

        if (it_time+1 == time.end())
            next_time = tmax;
        else
            next_time = *(it_time+1);

        for(auto it_edge = it_edges->begin();
                 it_edge != it_edges->end();
                 ++it_edge)
        {
            if(verbose)
                cout << "evaluating new edge" << endl;

            pair < size_t, size_t > &edge = *it_edge;

            if (edge.first > edge.second)
                swap(edge.first,edge.second);

            size_t edge_int = get_edge_integer(N,edge,hash_to_int);

            if (edge_int == these_edge_trajectories.size())
            {
                edge_trajectory_entry this_entry;
                this_entry.edge = edge;
                this_entry.time_pairs.push_back(make_pair(t,next_time));
                these_edge_trajectories.push_back(this_entry);
            }
            else
            {
                edge_trajectory_entry &this_entry = these_edge_trajectories[edge_int];
                if (this_entry.time_pairs.back().second == t)
                    this_entry.time_pairs.back().second = next_time;
                else
                    this_entry.time_pairs.push_back(make_pair(t,next_time));
            }

            if(verbose)
                cout << "this edge = ( " << edge.first << ", " << edge.second << " )" << endl;

        }

        if (return_edge_similarities)
        {
            for(auto it_edge = it_edges->begin();
                     it_edge != it_edges->end();
                     ++it_edge)
            {
                if(verbose)
                    cout << "evaluating new edge" << endl;

                pair < size_t, size_t > &edge = *it_edge;

                if (edge.first > edge.second)
                    swap(edge.first,edge.second);

                size_t edge_int = get_edge_integer(N,edge,hash_to_int);

                auto it_next_edge = it_edge+1;
                for(; it_next_edge != it_edges->end(); ++it_next_edge) 
                {

                    pair < size_t, size_t > & next_edge = *it_next_edge;

                    if(verbose)
                        cout << "    next edge = ( " << next_edge.first << ", " << next_edge.second << " )" << endl;

                    if (next_edge.first > next_edge.second)
                        swap(next_edge.first, next_edge.second);

                    size_t &i = edge.first;
                    size_t &j = edge.second;
                    size_t &u = next_edge.first;
                    size_t &v = next_edge.second;

                    if ( (i==u) or (i==v) or (j==u) or (j==v) )
                    {

                        if (verbose)
                            cout << "    now integrating edge pair" << endl;

                        size_t edge_int_2 = get_edge_integer(N,next_edge,hash_to_int);

                        if (verbose)
                            cout << "    eint1 = " << edge_int << ", eint2 = " << edge_int_2 << endl;

                        similarities[ 
                                        get_sorted_pair( edge_int, edge_int_2 ) 
                                    ].value += (next_time - t);

                        if (verbose)
                            cout << "    new_similarity = " << 
                                similarities[ 
                                                get_sorted_pair( edge_int, edge_int_2 ) 
                                            ].value << endl;

                    }
                }
            }
        }

        it_edges++;
        it_time++;
    }

    edge_trajectories traj;

    vector < tuple < size_t, size_t, double > > sim;

    for(auto const &entry: similarities)
    {
        sim.push_back( make_tuple( entry.first.first, entry.first.second, entry.second.value ));
    }
    
    traj.trajectories = these_edge_trajectories;
    traj.edge_similarities = sim;

    traj.N = list_of_edge_lists.N;
    traj.t0 = list_of_edge_lists.t[0];
    traj.tmax = list_of_edge_lists.tmax;
    traj.time_unit = list_of_edge_lists.time_unit;
    traj.notes = list_of_edge_lists.notes;
    traj.int_to_node = list_of_edge_lists.int_to_node;

    return traj;
    
}

edge_changes
    edge_trajectories_to_edge_changes(
                edge_trajectories &traj
            )
{
    // this map keeps a sorted list of edge_in/edge_out events.
    map < double, event_map_entry > events;

    // copy other properties
    double t0 = traj.t0;
    double tmax = traj.tmax;
    size_t N = traj.N;

    // initial edge list
    vector < pair < size_t, size_t > > edges_initial;

    // iterate through trajectories
    for( auto const &entry: traj.trajectories )
    {
        // fetch edge
        pair < size_t, size_t > edge = entry.edge;

        // iterate through time intervals of this edge
        for (auto const &interval : entry.time_pairs)
        {
            double const &ta = interval.first;
            double const &tb = interval.second;

            if (ta == t0)
                edges_initial.push_back(edge);
            else
                events[ta].edges_in.push_back(edge);

            if (tb != tmax)
                events[tb].edges_out.push_back(edge);
        }
    }

    vector < double > t;
    vector < vector < pair < size_t, size_t > > > edges_in;
    vector < vector < pair < size_t, size_t > > > edges_out;

    // iterate through sorted list
    for( auto & event: events )
    {
        t.push_back(event.first);
        edges_in.push_back(event.second.edges_in);
        edges_out.push_back(event.second.edges_out);
    }

    edge_changes ec;

    ec.N = N;
    ec.t0 = t0;
    ec.t = t;
    ec.tmax = tmax;
    ec.edges_initial = edges_initial;
    ec.edges_in = edges_in;
    ec.edges_out = edges_out;
    ec.time_unit = traj.time_unit;
    ec.notes = traj.notes;
    ec.int_to_node = traj.int_to_node;

    return ec;
}
