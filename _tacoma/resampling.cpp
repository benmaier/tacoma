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
#include "resampling.h"

using namespace std;

edge_lists
     sample_from_edge_lists(
             edge_lists &list_of_edge_lists,
             double dt,
             size_t N_time_steps,
             const bool sample_aggregates,
             const bool verbose
        )
{

    // get references to edge_list and time
    vector < vector < pair < size_t, size_t > > > & all_edges = list_of_edge_lists.edges;
    vector < double > time = list_of_edge_lists.t; // this is a copy! because we're pushing tmax at the end
    double t0 = time[0];
    double tmax = list_of_edge_lists.tmax;

    if (tmax < time.back())
        throw domain_error("The value tmax is smaller than the last value in the time list.");

    time.push_back(tmax);
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

    vector < vector < pair < size_t, size_t > > > new_list_of_edge_lists;
    vector < double > new_time;

    if( verbose)
    {
        cout << "starting resampling process with dt = " << dt << " and N_time_steps = " << N_time_steps << endl;
        cout << "sample_aggregates = " << sample_aggregates << endl;
    }

    for(size_t t_i = 0; t_i < N_time_steps; t_i++)
    {
        double this_new_time = dt * t_i + t0;

        if (sample_aggregates)
        {
            graph_from_edgelist(G,*it_edge_lists);
        }

        // while the next (upcoming) time value in the original data
        // is smaller than the new demanded time
        if (verbose)
        {
            cout << "current demanded new time is " << this_new_time << endl;
            cout << "current original time is " << (*it_time) << endl;
            cout << "upcoming original time is " << (*(it_time+1)) << endl;
        }
        while (it_time != time.end() and ((*(it_time+1)) <= this_new_time))
        {
            if (verbose)
                cout << "advancing original time" << endl;
            it_time++;
            it_edge_lists++;

            if (verbose)
                cout << "new original time is " << *it_time << endl;
            if (sample_aggregates)
            {
                for(auto const &edge: *it_edge_lists)
                {
                    G[ edge.first ].insert( edge.second );
                    G[ edge.second ].insert( edge.first );
                }
            }

        }

        if (not sample_aggregates)
            new_list_of_edge_lists.push_back(*it_edge_lists);
        else
        {
            vector < pair < size_t, size_t > > this_new_edge_list;
            edgelist_from_graph(this_new_edge_list, G);
            new_list_of_edge_lists.push_back(this_new_edge_list);
        }

        new_time.push_back(this_new_time);
    }

    edge_lists result;
    result.tmax = tmax;
    result.N = N;
    result.edges = new_list_of_edge_lists;
    result.t = new_time;

    return result;
}


edge_lists
     sample_from_edge_changes(
             edge_changes &list_of_edge_changes,
             double dt,
             size_t N_time_steps,
             const bool sample_aggregates,
             const bool verbose
             )
{
    // get references to edge_list and time
    vector < vector < pair < size_t, size_t > > > & all_edges_in = list_of_edge_changes.edges_in;
    vector < vector < pair < size_t, size_t > > > & all_edges_out = list_of_edge_changes.edges_out;
    vector < double > time = list_of_edge_changes.t; //this is a copy
    // set initial and final time
    double t0 = list_of_edge_changes.t0;
    double tmax = list_of_edge_changes.tmax;
    time.push_back(tmax);

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

    vector < vector < pair < size_t, size_t > > > new_list_of_edge_lists;
    vector < double > new_time;
    new_time.push_back(t0);
    new_list_of_edge_lists.push_back(list_of_edge_changes.edges_initial);

    if (tmax < time.back())
        throw domain_error("The value tmax is smaller than the last value in the time list.");

    if( verbose)
    {
        cout << "starting resampling process with dt = " << dt << " and N_time_steps = " << N_time_steps << endl;
        cout << "sample_aggregates = " << sample_aggregates << endl;
    }

    // start at new_time = t0 + dt because t0 was already taken care of
    for(size_t t_i = 1; t_i < N_time_steps; t_i++)
    {
        double this_new_time = dt * t_i + t0;

        vector < set < size_t > > G_aggregate;

        if (sample_aggregates)
        {
            G_aggregate = G;
        }

        if (verbose)
        {
            cout << "current demanded new time is " << this_new_time << endl;
            cout << "current original time is " << (*it_time) << endl;
            cout << "upcoming original time is " << (*(it_time+1)) << endl;
        }

        // while the next (upcoming) time value in the original data
        // is smaller than the new demanded time
        while (it_time != time.end() and ((*(it_time)) <= this_new_time))
        {
            if (verbose)
                cout << "advancing original time" << endl;
            if (verbose)
                cout << "new original time is " << *it_time << endl;
            for(auto const &edge: *it_edges_in)
            {

                size_t const & i = edge.first;
                size_t const & j = edge.second;
                if (sample_aggregates)
                {
                    G_aggregate[ i ].insert( j );
                    G_aggregate[ j ].insert( i );
                }
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
            
            it_time++;
            it_edges_in++;
            it_edges_out++;


        }

        vector < pair < size_t, size_t > > this_new_edge_list;

        if (sample_aggregates)
            edgelist_from_graph(this_new_edge_list,G_aggregate);
        else
            edgelist_from_graph(this_new_edge_list,G);

        new_list_of_edge_lists.push_back(this_new_edge_list);
        new_time.push_back(this_new_time);
    }

    edge_lists result;
    result.tmax = tmax;
    result.N = N;
    result.edges = new_list_of_edge_lists;
    result.t = new_time;

    return result;
}

edge_lists
     bin_from_edge_lists(
             edge_lists &list_of_edge_lists,
             double dt,
             size_t N_time_steps,
             const bool verbose
        )
{

    // get references to edge_list and time
    vector < vector < pair < size_t, size_t > > > & all_edges = list_of_edge_lists.edges;
    vector < double > & time = list_of_edge_lists.t;
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

    vector < vector < pair < size_t, size_t > > > new_list_of_edge_lists;
    vector < double > new_time;

    if( verbose)
    {
        cout << "starting binning process with dt = " << dt << " and N_time_steps = " << N_time_steps << endl;
    }


    for(size_t t_i = 0; t_i < N_time_steps; t_i++)
    {
        double this_new_time = dt * t_i + t0;
        double this_next_time = dt * (t_i+1) + t0;

        if (verbose)
            cout << "loading edges from original time " << *it_time << endl;

        graph_from_edgelist(G,*it_edge_lists);

        // while the next (upcoming) time value in the original data
        // is smaller than the new demanded time
        if (verbose)
        {
            cout << "loaded edges [ ";
            for(auto const &edge: *it_edge_lists)
                cout << "( " << edge.first << " " << edge.second << " ) ";
            cout << "]" << endl;
            cout << "current demanded new time is " << this_new_time << endl;
            cout << "current original time is " << (*it_time) << endl;
            cout << "upcoming demanded new time is " << this_next_time << endl;
            cout << "upcoming original time is " << (*(it_time+1)) << endl;
        }

        while ((it_time+1) != time.end() and (*(it_time+1) <= this_next_time))
        {
            it_time++;
            it_edge_lists++;

            if(verbose)
            {
                cout << "    adding edges [ ";
                for(auto const &edge: *it_edge_lists)
                    cout << "( " << edge.first << " " << edge.second << " ) ";
                cout << "]" << endl;
            }

            if (*(it_time) < this_next_time)
            {
                for(auto const &edge: *it_edge_lists)
                {
                    G[ edge.first ].insert( edge.second );
                    G[ edge.second ].insert( edge.first );
                }
            }

            if (verbose)
            {
                cout << "    advanced original time to " << *it_time << endl;
                cout << "    it_edge_lists == all_edges.end() = " << (it_edge_lists == all_edges.end()) << endl;
            }
        }

        
        if (verbose)
            cout << "left loop" << endl;               

        vector < pair < size_t, size_t > > this_new_edge_list;
        edgelist_from_graph(this_new_edge_list, G);
        new_list_of_edge_lists.push_back(this_new_edge_list);

        new_time.push_back(this_new_time);

    }

    edge_lists result;
    result.tmax = tmax;
    result.N = N;
    result.edges = new_list_of_edge_lists;
    result.t = new_time;

    return result;
}

edge_lists
     bin_from_edge_changes(
             edge_changes &list_of_edge_changes,
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

    vector < vector < pair < size_t, size_t > > > new_list_of_edge_lists;
    vector < double > new_time;

    if (tmax < time.back())
        throw domain_error("The value tmax is smaller than the last value in the time list.");

    if( verbose)
    {
        cout << "starting resampling process with dt = " << dt << " and N_time_steps = " << N_time_steps << endl;
    }

    // start at new_time = t0 + dt because t0 was already taken care of
    for(size_t t_i = 0; t_i < N_time_steps; t_i++)
    {
        double this_new_time = dt * t_i + t0;
        double this_next_time = dt * (t_i+1) + t0;

        vector < set < size_t > > G_aggregate = G;

        if (verbose)
        {
            cout << "current demanded new time is " << this_new_time << endl;
            cout << "current original time is " << (*it_time) << endl;
            cout << "upcoming original time is " << (*(it_time+1)) << endl;
        }

        // while the next (upcoming) time value in the original data
        // is smaller than the new demanded time
        while (it_time != time.end() and ((*(it_time)) <= this_next_time))
        {
            if (verbose)
                cout << "advancing original time" << endl;
            if (verbose)
                cout << "new original time is " << *it_time << endl;
            for(auto const &edge: *it_edges_in)
            {

                size_t const & i = edge.first;
                size_t const & j = edge.second;

                if (*it_time < this_next_time)
                {
                    G_aggregate[ i ].insert( j );
                    G_aggregate[ j ].insert( i );
                }
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
            
            it_time++;
            it_edges_in++;
            it_edges_out++;


        }

        vector < pair < size_t, size_t > > this_new_edge_list;

        edgelist_from_graph(this_new_edge_list,G_aggregate);

        new_list_of_edge_lists.push_back(this_new_edge_list);
        new_time.push_back(this_new_time);
    }

    edge_lists result;
    result.tmax = tmax;
    result.N = N;
    result.edges = new_list_of_edge_lists;
    result.t = new_time;

    return result;
}
