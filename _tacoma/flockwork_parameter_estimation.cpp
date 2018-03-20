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
#include "conversion.h"

using namespace std;


flockwork_args
     get_flockwork_P_args(
             edge_changes &list_of_edge_changes,
             double dt,
             size_t N_time_steps,
             double k_over_k_real_scaling,
             map < pair < size_t, size_t >, double > aggregated_network,
             const bool ensure_empty_network,
             const bool change_tmax_if_dt_does_not_fit,
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
            if (change_tmax_if_dt_does_not_fit)
            {
                N_time_steps = (size_t) ceil(_N_t);
                tmax = t0 + N_time_steps * dt;
            }
            else
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

    // get the observables for the binning process
    vector < double > new_time;
    vector < size_t > m_in;
    vector < size_t > m_out;
    vector < size_t > m;

    m.push_back( list_of_edge_changes.edges_initial.size() );
    new_time.push_back( t0 );

    if (tmax < time.back())
        throw domain_error("The value tmax is smaller than the last value in the time list.");

    if( verbose)
    {
        cout << "starting binning process with dt = " << dt << " and N_time_steps = " << N_time_steps << endl;
    }

    // start at new_time = t0 + dt because t0 was already taken care of
    for(size_t t_i = 1; t_i < N_time_steps; t_i++)
    {
        double this_new_time = dt * t_i + t0;

        // initialize the counts for this time bin
        m_in.push_back(0);
        m_out.push_back(0);

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
            {
                cout << "advancing original time" << endl;
                cout << "new original time is " << *it_time << endl;
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

            m_in.back() += it_edges_in->size();
            m_out.back() += it_edges_out->size();
            
            it_time++;
            it_edges_in++;
            it_edges_out++;
        }

        size_t num_edges = 0;

        for(auto const &neighbors: G)
            num_edges += neighbors.size();

        num_edges /= 2;

        m.push_back( num_edges );
        new_time.push_back( this_new_time );
    }

    // ============== add the last time frame for looping ================
    new_time.push_back( dt * N_time_steps + t0 );

    // get a sorted edge list
    vector < pair < size_t, size_t > > last_edge_list;
    edgelist_from_graph(last_edge_list,G);

    set < size_t > last_edge_integers = get_edge_integer_set(last_edge_list, N);
    set < size_t > edge_integers = get_edge_integer_set(list_of_edge_changes.edges_initial, N);

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

    m_in.push_back(incoming_edge_integers.size());
    m_out.push_back(outgoing_edge_integers.size());
    m.push_back( m.front() );

    // ========= now go through the observables and compute the parameters =============
    auto it_m = m.begin();
    auto it_m_in = m_in.begin();
    auto it_m_out = m_out.begin();
    it_time = new_time.begin();

    //initialize parameters
    vector < pair < double, double > > gamma;
    vector < double > P;

    while(it_m_in != m_in.end())
    {
        size_t & this_m = *it_m;
        size_t & this_m_in = *it_m_in;
        size_t & this_m_out = *it_m_out;

        if ( (this_m == 0) and ( this_m_out > 0) )
        {
            /* this happens, e.g. in the beginning for cumulated bins. originally,
               the network is a null graph, so no edges can go out.
               however, when binning, the first cumulative bin might contain
               multiple bins where there's already edges. so we trick a bit
               and fake and set the first m_in to m_in := m_in - m_out
               and the first m_out := 0
            */
            if (this_m_in < this_m_out)
                this_m_in = 0;
            else
                this_m_in -= this_m_out;

            this_m_out = 0;
        }

        double _P;
        double _g;

        double _m_out = (double) this_m_out;
        double _m_in = (double) this_m_in;
        double _m = (double) this_m;

        if (this_m_out > 0)
        {
            // standard case
            _g = 0.5 * _m_out / _m / dt;
            _P = min(
                     2.0 * _m / (N+2.0*_m) * _m_in / _m_out ,
                     1.0
                     );
        }
        else if ( (this_m_out == 0) and (this_m_in > 0) )
        {
            // secondary case if there's no outgoing edges
            double next_m = (double) *(it_m+1);
            _P = 2.0 * next_m / (N+2.0*next_m);
            _g = _m_in / ( _P * (N+2.0*_m) * dt);
        }
        else
        {
            // tertiary case
            _g = 0.0;
            _P = 0.0;
        }

        // if in the beginning of this bin and at the end of this bin
        // there's no 
        if ( ( ensure_empty_network ) and 
             ( this_m == 0 ) and 
             ( *(it_m+1) == 0 ) and
             ( _g == 0.0 ) and 
             ( _P == 0.0 )
           )
        {
            _g = 20.0 / dt / ((double) N);
        }

        gamma.push_back( make_pair( *it_time, _g * k_over_k_real_scaling) );
        P.push_back( _P / k_over_k_real_scaling );

        it_m++; 
        it_m_in++;
        it_m_out++; 
        it_time++;
    }

    flockwork_args fw_args;

    fw_args.N = N;
    fw_args.E = list_of_edge_changes.edges_initial;
    fw_args.P = P;
    fw_args.rewiring_rate = gamma;
    fw_args.tmax = tmax;
    fw_args.new_time = new_time;
    fw_args.m_in = m_in;
    fw_args.m_out = m_out;
    fw_args.m = m;

    if (aggregated_network.size() > 0)
    {
        vector < pair < vector < size_t >, vector < double > > > neighbor_affinity(N);

        for (auto const & entry: aggregated_network)
        {
            size_t i = entry.first.first;
            size_t j = entry.first.second;
            double val = entry.second;

            neighbor_affinity[i].first.push_back( j );
            neighbor_affinity[i].second.push_back( val );

            neighbor_affinity[j].first.push_back( i );
            neighbor_affinity[j].second.push_back( val );
        }

        fw_args.neighbor_affinity = neighbor_affinity;
    }
    else
    {
        vector < pair < vector < size_t >, vector < double > > > neighbor_affinity;
        fw_args.neighbor_affinity = neighbor_affinity;
    }

    return fw_args;
}

