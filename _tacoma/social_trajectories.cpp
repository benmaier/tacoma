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
#include "conversion.h"

using namespace std;

size_t hash_edge(const pair < size_t, size_t > &p, const size_t &N)
{
    return N * p.first + p.second;
}

size_t get_edge_integer(
            const size_t &N,
            pair < size_t, size_t > const &edge,
            map < size_t, size_t > &hash_to_int
        )
{
    size_t N_hashes = hash_to_int.size();

    size_t this_hash = hash_edge(edge,N);

    auto current_hash_pair = hash_to_int.find(this_hash);
    if (current_hash_pair != hash_to_int.end())
        return current_hash_pair->second;
    else
    {
        hash_to_int[this_hash] = N_hashes;
        return N_hashes;
    }
}

size_t get_group_integer(
            size_t const &N,
            set < size_t > const &component,
            map < size_t, size_t > &hash_to_int,
            size_t &this_hash
        )
{
    size_t N_hashes = hash_to_int.size();
    vector < bool > this_group(N, false);

    for( auto const &node: component)
    {
        this_group[node] = true;
    }

    hash < vector < bool > > hash_calculator;

    this_hash = hash_calculator(this_group);

    auto current_hash_pair = hash_to_int.find(this_hash);
    if (current_hash_pair != hash_to_int.end())
        return current_hash_pair->second;
    else
    {
        hash_to_int[this_hash] = N_hashes;
        return N_hashes;
    }
}

vector < set < size_t > > 
    binned_social_trajectory_from_edge_lists(
             edge_lists &list_of_edge_lists,
             size_t node,
             double dt,
             size_t N_time_steps,
             const bool verbose
        )
{

    // get references to edge_list and time
    vector < vector < pair < size_t, size_t > > > & all_edges = list_of_edge_lists.edges;
    vector < double > time = list_of_edge_lists.t;
    double t0 = time[0];
    double tmax = list_of_edge_lists.tmax;

    if (tmax < time.back())
        throw domain_error("The value tmax is smaller than the last value in the time list.");

    size_t N = list_of_edge_lists.N;

    // check w
    if ((dt > 0.0) and N_time_steps > 0)
        throw domain_error("please provide either positive dt or positive N_time_steps, not both positive");    
    else if ((dt == 0.0) and N_time_steps == 0)
        throw domain_error("please provide either positive dt or positive N_time_steps, not both zero");
    else if ((dt > 0.0) and N_time_steps == 0)
    {
        double _N_t = (tmax - time.front()) / dt;
        double intpart;

        // check if this yielded a nice round number
        // by checking the rest after the period
        if (modf(_N_t,&intpart) != 0.0)
            throw domain_error("dt does not nicely divide time interval (tmax - t0) in integer parts");
        else
            N_time_steps = (size_t) _N_t;
    }
    else if ((dt == 0.0) and N_time_steps > 0)
        dt = (tmax - time.front()) / N_time_steps;

    // create graph
    vector < set < size_t > > G(N);

    // create iteratos
    auto it_edge_lists = all_edges.begin();
    auto it_time = time.begin();

    vector < set < size_t > > trajectory(N_time_steps);
    map < size_t, size_t > hash_to_int;

    auto current_trajectory_bin = trajectory.begin();
    size_t old_bin_number = 0;
    size_t current_group_integer = 1;

    if( verbose)
    {
        cout << "starting resampling process with dt = " << dt << " and N_time_steps = " << N_time_steps << endl;
    }

    while (it_time != time.end())
    {
        double this_time = *it_time;
        size_t bin_number = (this_time-t0) / (tmax-t0) * N_time_steps;
        size_t bin_difference = bin_number - old_bin_number;

        if (verbose)
        {
            cout << "time of event = " << this_time << endl;
            cout << "bin_number = " << bin_number << endl;
            cout << "old_bin_number = " << old_bin_number << endl;
            cout << "bin_difference = " << bin_difference << endl;

        }

        for(size_t bin = 1; bin <= bin_difference; bin++)
        {
            current_trajectory_bin++;
            if (G[node].size() > 0 and (t0 + bin_number * dt > this_time) )
            {
                if (current_group_integer < hash_to_int.size() )
                    current_trajectory_bin->insert(current_group_integer);
            }
        }


        old_bin_number = bin_number;

        if (verbose)
        {
            cout << "advanced bins to bin " << old_bin_number << endl;
        }

        graph_from_edgelist(G,*it_edge_lists);
        
        if (verbose)
        {
            cout << "got next graph with neighbor list";
            for (auto const & neigh: G[node])
                cout << " " << neigh;
            cout << endl;
        }


        if (G[node].size() > 0)
        {
            set < size_t > component = get_component_of_node(node,G);
            size_t this_hash;
            current_group_integer = get_group_integer(N,
                                                      component,
                                                      hash_to_int,
                                                      this_hash
                                                     );
            if (verbose)
            {
                cout << "the current group integer is " << current_group_integer << endl;
            }

            current_trajectory_bin->insert(
                                            current_group_integer
                                          );
        }

        it_time++;
        it_edge_lists++;

    }

    return trajectory;
}

vector < set < size_t > >
    binned_social_trajectory_from_edge_changes(
             edge_changes &list_of_edge_changes,
             size_t node,
             double dt,
             size_t N_time_steps,
             const bool verbose
        )
{
    // get references to edge_list and time
    vector < vector < pair < size_t, size_t > > > & all_edges_in = list_of_edge_changes.edges_in;
    vector < vector < pair < size_t, size_t > > > & all_edges_out = list_of_edge_changes.edges_out;
    vector < double > & time = list_of_edge_changes.t;
    // set initial and final time
    double t0 = list_of_edge_changes.t0;
    double tmax = list_of_edge_changes.tmax;

    // create graph
    size_t N = list_of_edge_changes.N;
    vector < set < size_t > > G(N);
    graph_from_edgelist(G,list_of_edge_changes.edges_initial);

    // check w
    if ((dt > 0.0) and N_time_steps > 0)
        throw domain_error("please provide either positive dt or positive N_time_steps, not both positive");    
    else if ((dt == 0.0) and N_time_steps == 0)
        throw domain_error("please provide either positive dt or positive N_time_steps, not both zero");
    else if ((dt > 0.0) and N_time_steps == 0)
    {
        double _N_t = (tmax - t0) / dt;
        double intpart;

        // check if this yielded a nice round number
        // by checking the rest after the period
        if (modf(_N_t,&intpart) != 0.0)
            throw domain_error("dt does not nicely divide time interval (tmax - t0) in integer parts");
        else
            N_time_steps = (size_t) _N_t;
    }
    else if ((dt == 0.0) and N_time_steps > 0)
        dt = (tmax - t0) / N_time_steps;

    // create iteratos
    auto it_edges_in = all_edges_in.begin();
    auto it_edges_out = all_edges_out.begin();
    auto it_time = time.begin();

    vector < set < size_t > > trajectory(N_time_steps);
    map < size_t, size_t > hash_to_int;
    auto current_trajectory_bin = trajectory.begin();
    size_t old_bin_number = 0;


    size_t current_group_integer = 1;
    if (G[node].size() > 0)
    {
        set < size_t > component = get_component_of_node(node,G);
        size_t this_hash;
        current_group_integer = get_group_integer(N,
                                               component,
                                               hash_to_int,
                                               this_hash
                                               );
        current_trajectory_bin->insert(current_group_integer);
    }

    if (tmax < time.back())
        throw domain_error("The value tmax is smaller than the last value in the time list.");

    if( verbose)
    {
        cout << "starting resampling process with dt = " << dt << " and N_time_steps = " << N_time_steps << endl;
    }

    while (it_time != time.end())
    {
        double this_time = *it_time;
        size_t bin_number = (this_time-t0) / (tmax-t0) * N_time_steps;
        size_t bin_difference = bin_number - old_bin_number;

        if (verbose)
        {
            cout << "time of event = " << this_time << endl;
            cout << "bin_number = " << bin_number << endl;
            cout << "old_bin_number = " << old_bin_number << endl;
            cout << "bin_difference = " << bin_difference << endl;

        }

        for(size_t bin = 1; bin <= bin_difference; bin++)
        {
            if (verbose)
            {
                cout << "advancing to bin " << old_bin_number + bin << endl;
            }
            current_trajectory_bin++;
            if (G[node].size() > 0 and (t0 + bin_number * dt > this_time) )
            {
                if (current_group_integer < hash_to_int.size() )
                {
                    current_trajectory_bin->insert(current_group_integer);
                    if (verbose)
                    {
                        cout << "added group integer " << current_group_integer << endl;
                    }
                }
            }
        }

        old_bin_number = bin_number;

        if (verbose)
        {
            cout << "advanced bins to bin " << old_bin_number << endl;
        }

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

        if (verbose)
        {
            cout << "got next graph with neighbor list";
            for (auto const & neigh: G[node])
                cout << " " << neigh;
            cout << endl;
        }

        if (G[node].size() > 0)
        {
            if (verbose)
            {
                cout << "the current group integer is " << current_group_integer << endl;
            }

            set < size_t > component = get_component_of_node(node,G);
            size_t this_hash;
            current_group_integer = get_group_integer(N,
                                                      component,
                                                      hash_to_int,
                                                      this_hash
                                                     );
            current_trajectory_bin->insert( current_group_integer );
        }

        it_time++;
        it_edges_in++;
        it_edges_out++;
    }

    return trajectory;
}

vector < social_trajectory_entry >
     social_trajectory_from_edge_changes(
            edge_changes &list_of_edge_changes,
            size_t node,
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

    if (verbose)
    {
        cout << "last time in array = " << time.back() << endl;
        cout << "              tmax = " << tmax << endl;
        cout << " created Graph with degree sequence ";
        for( auto const &neighbors: G)
        {
            cout << neighbors.size() << " ";
        }
        cout << endl;

    }

    if (tmax < time.back())
        throw domain_error("The value tmax is smaller than the last time in the time list.");

    // create iterators
    auto it_edges_in = all_edges_in.begin();
    auto it_edges_out = all_edges_out.begin();
    auto it_time = time.begin();

    // container to track activities of nodes
    double last_time_active = t0;

    // vectors and sets dealing with the change of components
    set < size_t > old_component = get_component_of_node(node,G);
    set < size_t > new_component;
    map < size_t, size_t > hash_to_int;

    vector < social_trajectory_entry > social_trajectory;

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

        if (verbose)
        {
            cout << " created Graph with degree sequence ";
            for( auto const &neighbors: G)
            {
                cout << neighbors.size() << " ";
            }
            cout << endl;
        }

        new_component = get_component_of_node(node,G);

        if (verbose)
        {
            cout << "found component containing nodes ";
            for( auto const &neighbor: new_component)
            {
                cout << neighbor << " ";
            }
            cout << endl;
        }

        bool component_changed = (old_component != new_component);

        if ((old_component.size()>1) and (component_changed) and (t>t0))
        {
            size_t this_hash;
            size_t this_group_int = get_group_integer(N,old_component,hash_to_int,this_hash);

            if (this_group_int == social_trajectory.size())
            {
                social_trajectory_entry this_entry;
                vector < pair < double, double > > times;

                this_entry.size = old_component.size();
                this_entry.hash = this_hash;
                this_entry.time_pairs = times;

                social_trajectory.push_back(this_entry);
            }

            social_trajectory[this_group_int].time_pairs.push_back( make_pair(last_time_active,t) );
        }

        if (component_changed)
        {
            last_time_active = t;
        }

        // copy this graph and the components for comparisons in next time slice
        old_component = new_component;

        // if this is the last time step, then not only save the potentially old component
        // but also the new component

        // advance iterators
        it_edges_in++;
        it_edges_out++;
        it_time++;

        if (it_time == time.end() and old_component.size()>1)
        {
            t = tmax;
            size_t this_hash;
            size_t this_group_int = get_group_integer(N,old_component,hash_to_int,this_hash);

            if (this_group_int == social_trajectory.size())
            {
                social_trajectory_entry this_entry;
                vector < pair < double, double > > times;

                this_entry.size = old_component.size();
                this_entry.hash = this_hash;
                this_entry.time_pairs = times;

                social_trajectory.push_back(this_entry);
            }

            social_trajectory[this_group_int].time_pairs.push_back( make_pair(last_time_active,t) );
            last_time_active = t;
        }
    }

    return social_trajectory;
    
}

vector < social_trajectory_entry >
     social_trajectory_from_edge_lists(
            edge_lists &list_of_edge_lists,
            size_t node,
            const bool verbose
        )
{
    // get references to edge_list and time
    vector < vector < pair < size_t, size_t > > > & all_edges = list_of_edge_lists.edges;
    vector < double > & time = list_of_edge_lists.t;

    // create graph
    size_t N = list_of_edge_lists.N;
    vector < set < size_t > > G(N);

    // set initial and final time
    double t0 = time.front();
    double tmax = list_of_edge_lists.tmax;

    if (verbose)
    {
        cout << "last time in array = " << time.back() << endl;
        cout << "              tmax = " << tmax << endl;
        cout << " created Graph with degree sequence ";
        for( auto const &neighbors: G)
        {
            cout << neighbors.size() << " ";
        }
        cout << endl;

    }

    if (tmax < time.back())
        throw domain_error("The value tmax is smaller than the last time in the time list.");

    // create iterators
    auto it_edges = all_edges.begin();
    auto it_time = time.begin();

    // container to track activities of nodes
    double last_time_active = t0;

    // vectors and sets dealing with the change of components
    set < size_t > old_component = get_component_of_node(node,G);
    set < size_t > new_component;
    map < size_t, size_t > hash_to_int;

    vector < social_trajectory_entry > social_trajectory;

    // iterate through all changes
    while (it_edges != all_edges.end() and it_time != time.end() )
    {
        // get time of current state
        double t = *it_time;

        // get a sorted edge list
        vector < pair < size_t, size_t > > & these_edges = (*it_edges);
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

        new_component = get_component_of_node(node,G);

        if (verbose)
        {
            cout << "found component containing nodes ";
            for( auto const &neighbor: new_component)
            {
                cout << neighbor << " ";
            }
            cout << endl;
        }

        bool component_changed = (old_component != new_component);

        if ((old_component.size()>1) and (component_changed) and (t>t0))
        {
            size_t this_hash;
            size_t this_group_int = get_group_integer(N,old_component,hash_to_int,this_hash);

            if (this_group_int == social_trajectory.size())
            {
                social_trajectory_entry this_entry;
                vector < pair < double, double > > times;

                this_entry.size = old_component.size();
                this_entry.hash = this_hash;
                this_entry.time_pairs = times;

                social_trajectory.push_back(this_entry);
            }

            social_trajectory[this_group_int].time_pairs.push_back( make_pair(last_time_active,t) );
        }

        if (component_changed)
        {
            last_time_active = t;
        }

        // copy this graph and the components for comparisons in next time slice
        old_component = new_component;

        // if this is the last time step, then not only save the potentially old component
        // but also the new component

        // advance iterators
        it_edges++;
        it_time++;

        if (it_time == time.end() and old_component.size()>1)
        {
            t = tmax;
            size_t this_hash;
            size_t this_group_int = get_group_integer(N,old_component,hash_to_int,this_hash);

            if (this_group_int == social_trajectory.size())
            {
                social_trajectory_entry this_entry;
                vector < pair < double, double > > times;

                this_entry.size = old_component.size();
                this_entry.hash = this_hash;
                this_entry.time_pairs = times;

                social_trajectory.push_back(this_entry);
            }

            social_trajectory[this_group_int].time_pairs.push_back( make_pair(last_time_active,t) );
            last_time_active = t;
        }
    }

    return social_trajectory;
    
}

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

    return traj;
    
}

