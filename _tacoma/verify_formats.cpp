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
#include "social_trajectories.h"
#include "verify_formats.h"

using namespace std;

//vector < edge_trajectory_entry >
void
    verify_edge_changes(
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

    // set initial and final time
    double t0 = list_of_edge_changes.t0;
    double tmax = list_of_edge_changes.tmax;

    if (verbose)
    {
        cout << "last time in array = " << time.back() << endl;
        cout << "              tmax = " << tmax << endl;
    }

    if (tmax <= time.back())
        //throw domain_error(
        cout << "The value tmax is smaller or equal than the last time in the time list."
            //);
            << endl;

    if (all_edges_in.size() != time.size() or all_edges_out.size() != time.size() or all_edges_in.size() != all_edges_out.size())
        cout << "size of `edges_in` = " + to_string(all_edges_in.size()) + " != size of `t` = " + to_string(time.size()) << " != size of `edges_out` = " << to_string(all_edges_out.size()) << endl;

    // create iterators
    auto it_edges_in = all_edges_in.begin();
    auto it_edges_out = all_edges_out.begin();
    auto it_time = time.begin();


    // vectors and sets dealing with the change of edges
    map < size_t, size_t > hash_to_int;

    vector < edge_trajectory_entry > edge_trajectories;

    for(auto &edge: list_of_edge_changes.edges_initial)
    {
        if (edge.first > edge.second)
            swap(edge.first,edge.second);

        if (edge.first == edge.second)
            //throw domain_error(
                    cout << 
                               "The edge ("+to_string(edge.first)+", "
                               +to_string(edge.second)+") in the initial edge list "
                               +"is a self-loop."
                               //);
                        << endl;

        if ((edge.first >= N) or (edge.second >= N))
            //throw domain_error(
                    cout <<
                               "The edge ("+to_string(edge.first)+", "
                               +to_string(edge.second)+") in the initial edge list "
                               +"contains a node larger or equal to the number of nodes N = "+to_string(N)
                               //);
                        << endl;

        size_t edge_int = get_edge_integer(N,edge,hash_to_int);
        edge_trajectory_entry this_entry;
        if (edge_int == edge_trajectories.size())
        {
            this_entry.edge = edge;
            this_entry.last_time_active = t0;
            this_entry.is_active = true;
        }
        else
            //throw domain_error(
                    cout << "Doubling of edges in edge_changes.edges_initial"
                    //);
                        << endl;

        edge_trajectories.push_back(this_entry);
    }

    size_t t_count = 0;
    double last_time = t0;

    // iterate through all changes
    while (it_edges_in != all_edges_in.end() and it_time != time.end() and it_edges_out != all_edges_out.end() )
    {
        // get time of current state
        double t = *it_time;
        if (t <= last_time)
            //throw domain_error(
             cout <<
                                "this time value t = " + to_string(t) + " at it = " + to_string(t_count)
                               +" is smaller or equal than"
                               +" the last time value t_last = " + to_string(last_time)
                               //);
                        << endl;

        // get a sorted edge list
        vector < pair < size_t, size_t > > & these_edges_in = (*it_edges_in);
        vector < pair < size_t, size_t > > & these_edges_out = (*it_edges_out);

        set < size_t > edge_int_in;
        set < size_t > edge_int_out;

        for(auto &edge: these_edges_out)
        {
            if (edge.first > edge.second)
                swap(edge.first,edge.second);

            if (edge.first == edge.second)
                //throw domain_error(
                        cout << "The edge ("+to_string(edge.first)+", "
                                   +to_string(edge.second)+") in the edges_out list at time it = "+to_string(t_count)
                                   +" is a self-loop."
                                   //);
                        << endl;

            if ((edge.first >= N) or (edge.second >= N))
                //throw domain_error(
                        cout << 
                                  "The edge ("+to_string(edge.first)+", "
                                   +to_string(edge.second)+") in the edges_out list at time it = "+to_string(t_count)
                                   +" contains a node larger or equal to the number of nodes N = "+to_string(N)
                                   //);
                            << endl;
            edge_int_out.insert( hash_edge(edge, N) );
        }

        for(auto &edge: these_edges_in)
        {
            if (edge.first > edge.second)
                swap(edge.first,edge.second);

            if (edge.first == edge.second)
                //throw domain_error(
                        cout <<
                                    "The edge ("+to_string(edge.first)+", "
                                   +to_string(edge.second)+") in the edges_in list at time it = "+to_string(t_count)
                                   +" is a self-loop."
                                   //);
                            << endl;

            if ((edge.first >= N) or (edge.second >= N))
                //throw domain_error(
                cout<<
                        "The edge ("+to_string(edge.first)+", "
                                   +to_string(edge.second)+") in the edges_in list at time it = "+to_string(t_count)
                                   +" contains a node larger or equal to the number of nodes N = "+to_string(N)
                                   //);
                            << endl;

            edge_int_in.insert( hash_edge(edge, N) );
        }

        vector < size_t > intersection;
        set_intersection(edge_int_in.begin(),edge_int_in.end(),
                         edge_int_out.begin(),edge_int_out.end(),
                         back_inserter(intersection));

        if (intersection.size() > 0)
        {
            for(auto const &e: intersection)
            {
                size_t j = e % N;
                size_t i = e / N;
                cout << "The edge ("<< i << ", " << j << ") is in both lists edges_in "
                     << "and edges_out at time it = " << t_count << endl;
            }
            //throw domain_error("Equal edges in both lists edges_in and edges_out. See output for details.");
        }


        for(auto &edge: these_edges_out)
        {
            size_t edge_int = get_edge_integer(N,edge,hash_to_int);

            if (edge_int == edge_trajectories.size())
                //throw domain_error(
                cout << 
                                   "The edge ("+to_string(edge.first)+", "
                                   +to_string(edge.second)+") in the edges_out "
                                   +"list at time it = "+to_string(t_count)
                                   +" did not exist before"
                                   //);
                            << endl;

            if ((edge.first >= N) or (edge.second >= N))
                //throw domain_error(
                cout<<
                        "The edge ("+to_string(edge.first)+", "
                                   +to_string(edge.second)+") in the edges_out list at time it = "+to_string(t_count)
                                   +" contains a node larger or equal to the number of nodes N = "+to_string(N)
                                   //);
                            << endl;

            edge_trajectory_entry &this_entry = edge_trajectories[edge_int];

            if (not this_entry.is_active)
                //throw domain_error(
                        cout << 
                                   "The edge ("+to_string(this_entry.edge.first)+", "
                                   +to_string(this_entry.edge.second)+") was supposed to be deactivated "
                                   +"at time it = "+to_string(t_count)
                                   +" but was not active before."
                                   //);
                            << endl;

            this_entry.time_pairs.push_back(make_pair(this_entry.last_time_active, t));
            this_entry.last_time_active = t0 - 1000.0;
            this_entry.is_active = false;
        }

        for(auto &edge: these_edges_in)
        {
            if (edge.first > edge.second)
                swap(edge.first,edge.second);

            size_t edge_int = get_edge_integer(N,edge,hash_to_int);

            if (edge_int == edge_trajectories.size())
            {
                edge_trajectory_entry this_entry;
                this_entry.edge = edge;
                this_entry.last_time_active = t;
                this_entry.is_active = true;
                edge_trajectories.push_back(this_entry);
            }
            else
            {
                edge_trajectory_entry & this_entry = edge_trajectories[edge_int];

                if (this_entry.is_active)
                    //throw domain_error(
                            cout << 
                                      "The edge ("+to_string(this_entry.edge.first)+", "
                                       +to_string(this_entry.edge.second)+") was supposed to be activated "
                                       +"at time it = "+to_string(t_count)
                                       +" but was already activated before at time t = "
                                       +to_string(this_entry.last_time_active)
                                       //);
                            << endl;

                this_entry.last_time_active = t;
                this_entry.is_active = true;
            }

        }

        it_edges_in++;
        it_edges_out++;
        it_time++;
        t_count++;
        last_time = t;
    }

    for(auto &this_entry: edge_trajectories)
    {
        if (this_entry.is_active)
            this_entry.time_pairs.push_back(make_pair(this_entry.last_time_active,tmax));
    }
}

//vector < edge_trajectory_entry >
void
        verify_edge_lists(
            edge_lists &list_of_edge_lists,
            const bool verbose
        )
{
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
        //throw domain_error(
                cout << "The value tmax is smaller than the last time in the time list."
                //); 
                << endl;

    if (all_edges.size() != time.size())
        cout << "size of `edges` = " + to_string(all_edges.size()) + " != size of `t` = " + to_string(time.size()) << endl;

    // create iterators
    auto it_edges = all_edges.begin();
    auto it_time = time.begin();

    size_t t_count = 0;
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

        if (next_time <= t)
        {
            cout << "the upcoming time t_{i+1} = " + to_string(next_time) + " is smaller or equal to the current time t_{i} = " + to_string(t) + " with i = " + to_string(t_count) << endl;
        }

        for(auto &edge: *it_edges)
        {
            if (edge.first > edge.second)
                swap(edge.first,edge.second);

            if (edge.first == edge.second)
                //throw domain_error(
                        cout << 
                                   "The edge ("+to_string(edge.first)+", "
                                   +to_string(edge.second)+") in the edge list at it = " + to_string(t_count)
                                   +" is a self-loop."
                                   //);
                             << endl;
            if ((edge.first >= N) or (edge.second >= N))
                //throw domain_error(
                cout<<
                        "The edge ("+to_string(edge.first)+", "
                                   +to_string(edge.second)+") in the edges list at time it = "+to_string(t_count)
                                   +" contains a node larger or equal to the number of nodes N = "+to_string(N)
                                   //);
                            << endl;
        }

        it_edges++;
        it_time++;
        t_count++;
    }
}
