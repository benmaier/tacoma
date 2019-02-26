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
#include "resampling.h"
#include "conversion.h"
#include "FW_P_varying.h"

using namespace std;


flockwork_args
     get_flockwork_P_args(
             edge_changes &list_of_edge_changes,
             double dt,
             size_t N_time_steps,
             double k_over_k_real_scaling,
             double gamma_scaling,
             double P_scaling,
             map < pair < size_t, size_t >, double > &aggregated_network,
             const bool ensure_empty_network,
             const bool adjust_last_bin_if_dt_does_not_fit,
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
            if (adjust_last_bin_if_dt_does_not_fit)
            {
                // if it did not, change the time of the network
                // to be larger than originally.
                // This will be changed back later when 
                // computing gamma and P
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

    bool last_bin_was_for_emptying = false;

    while(it_m_in != m_in.end())
    {
        size_t & this_m = *it_m;
        size_t & this_m_in = *it_m_in;
        size_t & this_m_out = *it_m_out;

        // if tmax was shifted because dt did not fit, rescale the last
        // time dt to the appropriate dt
        if (it_m_in + 1 == m_in.end())
        {
            double this_time = *it_time;
            double next_time = list_of_edge_changes.tmax;
            dt = next_time - this_time;
        }


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
        double _m = ((double) this_m) / k_over_k_real_scaling;

        if (this_m_out > 0)
        {
            // standard case
            _g = 0.5 * _m_out / _m / dt * gamma_scaling;
            _P = min(
                     2.0 * _m / (N+2.0*_m) * _m_in / _m_out * P_scaling,
                     1.0
                     );
        }
        else if ( (this_m_out == 0) and (this_m_in > 0) )
        {
            // secondary case if there's no outgoing edges
            double next_m = ((double) *(it_m+1)) / k_over_k_real_scaling;
            _P = 2.0 * next_m / (N+2.0*next_m) * P_scaling;
            _g = _m_in / ( _P * (N+2.0*_m) * dt) * gamma_scaling;
        }
        else
        {
            // tertiary case
            _g = 0.0;
            _P = 0.0;
        }

        // if in the beginning of this bin and at the end of this bin
        // there's no edges and we want this to be ensured
        if ( ( ensure_empty_network ) and 
             ( this_m == 0 ) and 
             ( *(it_m+1) == 0 ) and
             ( _g == 0.0 ) and 
             ( _P == 0.0 )
           )
        {
            // we want to pick every node in the network at least once, s.t.
            // all nodes have zero edges in the end of the time bin.
            // This is the classic problem of "how often do I have to pick
            // in order to have picked every number at least once" with solution
            // <T> = N log(N)
            // Since we want <Events([t,t+dt])> = N * log(N), we need
            //  gamma_per_node = (N log(N)) / (N / dt) such that
            //  \int_t^{t+dt} dt' gamma_per_node * N = N * log(N)
            //
            //  Also, instead of doing this for two time bins, we double the rate
            //  and just do it for one time bin
            if ( not last_bin_was_for_emptying )
                _g = 1 * log(N) / dt;

            last_bin_was_for_emptying = true;
        }
        else
        {
            last_bin_was_for_emptying = false;
        }

        gamma.push_back( make_pair( *it_time, _g ) );
        P.push_back( _P );

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
    fw_args.tmax = list_of_edge_changes.tmax;
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

double k_simulated_over_k_real(
                vector < pair < double, double > > &k_original, 
                vector < pair < double, double > > &k_simulated
             )
{

    if (k_simulated.size() != k_original.size())
        throw length_error("Original network and simulated network differ in number of time points.");

    double sum = 0.0;
    size_t N_meas = 0;

    auto sim = k_simulated.begin();
    auto orig = k_original.begin();
    for(; sim != k_simulated.end(); ++sim)
    {
        double &k_s = sim->second;
        double &k_o = orig->second;

        if ((k_o > 0.0) and (k_s > 0.0))
        {
            sum += k_s / k_o;
            N_meas += 1;
        }
        //else if ((k_o == 0.0) and (k_s == 0.0))
        //{
        //    sum += 1.0;
        //    N_meas += 1;
        //}


        ++orig;
    }

    return sum / (double) N_meas;
}

double k_RMSE(
             edge_lists &original_binned,
             edge_lists &simulated_binned
             )
{
    vector < pair < double, double > > k_simulated = mean_degree_from_edge_lists(simulated_binned);
    vector < pair < double, double > > k_original = mean_degree_from_edge_lists(original_binned);

    if (k_simulated.size() != k_original.size())
        throw length_error("Original network and simulated network differ in number of time points.");

    double sum = 0.0;
    size_t N_meas = 0;

    auto sim = k_simulated.begin();
    auto orig = k_original.begin();
    for(; sim != k_simulated.end(); ++sim)
    {
        double &k_s = sim->second;
        double &k_o = orig->second;
        sum += pow( k_s - k_o, 2);

        ++orig;
        ++N_meas;
    }

    return sqrt(sum / (double) N_meas);
}

double estimate_k_scaling_gradient_descent(
             edge_changes &original_edge_changes,
             double dt_for_inference,
             double dt_for_binning,
             size_t measurements_per_configuration,
             double learning_rate,
             double relative_error,
             size_t N_eval_max,
             bool   verbose
        )
{
    edge_lists original_binned = bin_from_edge_changes(original_edge_changes, dt_for_binning);
    vector < pair < double, double > > k_original = mean_degree_from_edge_lists(original_binned);

    double this_k_over_k_real_scaling = 1.0;
    double last_k_over_k_real_scaling = 0.0;

    size_t N_eval = 0;

    size_t update_every = 1;
    map < pair < size_t, size_t >, double > dummy_network;

    while ( 
            ( abs(last_k_over_k_real_scaling - this_k_over_k_real_scaling) / this_k_over_k_real_scaling > relative_error )
           and
            ( N_eval < N_eval_max )
          )
    {

        double mean = 0.0;
        vector < double > vals;
        vector < pair < double, double > > new_k( k_original.size() );

        for(unsigned int measurement = 0; measurement < measurements_per_configuration; ++measurement)
        {
            flockwork_args these_args = get_flockwork_P_args(
                                            original_edge_changes,
                                            dt_for_inference,
                                            0,
                                            this_k_over_k_real_scaling,
                                            1.0,
                                            1.0,
                                            dummy_network,
                                            true,
                                            true,
                                            false
                                            );
            edge_changes fw = flockwork_P_varying_rates(
                                these_args.E,
                                these_args.N,
                                these_args.P,
                                these_args.tmax,
                                these_args.rewiring_rate,
                                these_args.tmax
                              );

            edge_lists fw_binned = bin_from_edge_changes(fw, dt_for_binning);

            vector < pair < double, double > > this_k = mean_degree_from_edge_lists(fw_binned);
            auto this_t_k = this_k.begin();
            for(auto &t_k: new_k)
            {
                t_k.second += this_t_k->second / (double) measurements_per_configuration;
                ++this_t_k;
            }

        }

        mean = k_simulated_over_k_real(k_original, new_k);

        last_k_over_k_real_scaling = this_k_over_k_real_scaling;

        double inverse_scaling = 1.0 / this_k_over_k_real_scaling - learning_rate * (1.0-1.0/mean);
        this_k_over_k_real_scaling = 1.0 / inverse_scaling;

        if (verbose and N_eval % update_every == 0)
        {
            cout << "== generation: " << N_eval+1 << " ==" << endl;
            cout << "  k_simulated_over_k_real = " << mean << endl;
            cout << "  k_scaling = " << this_k_over_k_real_scaling << endl;
            cout << "  scaling_err = " << abs(last_k_over_k_real_scaling - this_k_over_k_real_scaling) / this_k_over_k_real_scaling << endl;
        }

        ++N_eval;
    }

    return this_k_over_k_real_scaling;
}

double estimate_k_scaling_gradient_descent_RMSE(
             edge_changes &original_edge_changes,
             double dt_for_inference,
             double dt_for_binning,
             size_t measurements_per_configuration,
             double learning_rate,
             double relative_error,
             size_t N_eval_max,
             bool   verbose
        )
{
    edge_lists original_binned = bin_from_edge_changes(original_edge_changes, dt_for_binning);

    double this_k_over_k_real_scaling = 1.0;
    double last_k_over_k_real_scaling = 0.0;

    size_t N_eval = 0;

    size_t update_every = 1;
    map < pair < size_t, size_t >, double > dummy_network;

    double last_RMSE;
    double d_scaling = 0.1;

    double mean = 0.0;
    vector < double > vals;

    for(unsigned int measurement = 0; measurement < measurements_per_configuration; ++measurement)
    {
        flockwork_args these_args = get_flockwork_P_args(
                                        original_edge_changes,
                                        dt_for_inference,
                                        0,
                                        this_k_over_k_real_scaling,
                                        1.0,
                                        1.0,
                                        dummy_network,
                                        true,
                                        true,
                                        false
                                        );
        edge_changes fw = flockwork_P_varying_rates(
                            these_args.E,
                            these_args.N,
                            these_args.P,
                            these_args.tmax,
                            these_args.rewiring_rate,
                            these_args.tmax
                          );

        edge_lists fw_binned = bin_from_edge_changes(fw, dt_for_binning);

        double this_value;

        this_value = k_RMSE(original_binned, fw_binned);

        mean += this_value;
        vals.push_back( this_value );
    }
    mean /= (double) measurements_per_configuration;
    last_RMSE = mean;

    this_k_over_k_real_scaling = 1.0 / (1.0 / this_k_over_k_real_scaling - d_scaling);

    while ( 
            ( abs(last_k_over_k_real_scaling - this_k_over_k_real_scaling) / this_k_over_k_real_scaling > relative_error )
           and
            ( N_eval < N_eval_max )
          )
    {

        double mean = 0.0;
        vector < double > vals;

        for(unsigned int measurement = 0; measurement < measurements_per_configuration; ++measurement)
        {
            flockwork_args these_args = get_flockwork_P_args(
                                            original_edge_changes,
                                            dt_for_inference,
                                            0,
                                            this_k_over_k_real_scaling,
                                            1.0,
                                            1.0,
                                            dummy_network,
                                            true,
                                            true,
                                            false
                                            );
            edge_changes fw = flockwork_P_varying_rates(
                                these_args.E,
                                these_args.N,
                                these_args.P,
                                these_args.tmax,
                                these_args.rewiring_rate,
                                these_args.tmax
                              );

            edge_lists fw_binned = bin_from_edge_changes(fw, dt_for_binning);

            double this_value;

            this_value = k_RMSE(original_binned, fw_binned);

            mean += this_value;
            vals.push_back( this_value );
        }
        mean /= (double) measurements_per_configuration;

        double err = 0.0;
        for(auto val = vals.begin(); val != vals.end(); ++val)
            err += pow(mean - *val,2);
        err /= (double) measurements_per_configuration*(measurements_per_configuration-1);
        err = sqrt(err);

        last_k_over_k_real_scaling = this_k_over_k_real_scaling;

        double inverse_scaling = 1.0 / this_k_over_k_real_scaling - d_scaling;
        this_k_over_k_real_scaling = 1.0 / inverse_scaling;

        if (verbose and N_eval % update_every == 0)
        {
            cout << "== generation: " << N_eval+1 << " ==" << endl;
            
            cout << "  k_RMSE = " << mean << " +/- " << err << endl;

            cout << "  k_scaling = " << this_k_over_k_real_scaling << endl;
            cout << "  scaling_err = " << abs(last_k_over_k_real_scaling - this_k_over_k_real_scaling) / this_k_over_k_real_scaling << endl;
        }

        ++N_eval;
    }

    return this_k_over_k_real_scaling;
}

pair < vector < double >, vector < double > >
    get_node_gamma_and_P(
                edge_changes &ec,
                vector < pair < double, double > > &gamma,
                vector < double > &P,
                bool use_event_rate_method
                )
{
    // get references to edge_list and time
    vector < vector < pair < size_t, size_t > > > & all_edges_in = ec.edges_in;
    vector < vector < pair < size_t, size_t > > > & all_edges_out = ec.edges_out;
    vector < double > &ec_time = ec.t;

    // get vectors from gamma and P
    vector < double > ab_time;
    vector < double > g;

    for(auto &g_entry: gamma)
    {
        ab_time.push_back(g_entry.first);
        g.push_back(g_entry.second);
    }

    // set initial and final time
    double tmax = ec.tmax;

    // graph observables
    size_t N = ec.N;
    double Nf(N);
    size_t number_of_edges = ec.edges_initial.size();

    // intialize the node properties
    vector < size_t > k_node(N);
    vector < size_t > M_out(N);
    vector < size_t > M_in(N);

    // intialize the integrals
    vector < double > I_out(N);
    vector < double > I_out_2(N);
    vector < double > I_in_1(N);
    vector < double > I_in_2(N);

    for(auto const &edge : ec.edges_initial)
    {
        ++k_node[edge.first];
        ++k_node[edge.second];
    }

    auto it_e_in = all_edges_in.begin();
    auto it_e_out = all_edges_out.begin();

    auto it_t_ab = ab_time.begin();
    auto it_t_ec = ec_time.begin();

    auto it_g = g.begin();
    auto it_P = P.begin();

    double t = ec.t0;

    double t_next_ab;
    double t_next_ec;
    double t_next;

    do
    {

        // check which is the next event time at which stuff changes
        if ((it_t_ab == ab_time.end()) or (it_t_ab+1 == ab_time.end()))
            t_next_ab = tmax;
        else
            t_next_ab = *(it_t_ab+1);

        if ((it_t_ec == ec_time.end()) or (it_t_ec+1 == ec_time.end()))
            t_next_ec = tmax;
        else
            t_next_ec = *(it_t_ec+1);

        if ( t_next_ec <= t_next_ab )
            t_next = t_next_ec;
        else
            t_next = t_next_ab;

        // compute the integrals up to the next change

        double dt = t_next - t;
        double mean_degree = ((double) number_of_edges) * 2.0 / Nf;

        auto i_k = k_node.begin();
        auto i_I_out = I_out.begin();
        auto i_I_out_2 = I_out_2.begin();
        auto i_I_in_1 = I_in_1.begin();
        auto i_I_in_2 = I_in_2.begin();

        double &_g = *it_g;
        double &_P = *it_P;
        double a = _g * _P;
        double b = _g * (1.0-_P);

        while(i_k != k_node.end())
        {
            double _k = (double) *i_k;

            if (use_event_rate_method)
            {
                (*i_I_out) += a * _k * dt;
                (*i_I_out_2) += b * _k * dt;
                (*i_I_in_1) += a * (Nf - _k - 1.0) * (_k + 1.0) / (Nf - 1.0) * dt;
                (*i_I_in_2) += a * (1.0 + mean_degree) * dt;
            } else {
                (*i_I_out) += _g * _k * dt;
                (*i_I_in_1) += _g * _P * (Nf - _k - 1.0) * (_k + 1.0) / (Nf - 1.0) * dt;
                (*i_I_in_2) += _g * _P * (1.0 + mean_degree) * dt;
            }

            ++i_k;
            ++i_I_out;
            ++i_I_out_2;
            ++i_I_in_1;
            ++i_I_in_2;
        }

        // advance time
        t = t_next;

        // change either the network or gamma and P or both

        if (t_next_ec <= t_next_ab)
        {
            // update observables for this advancement
            for(auto const &edge : *it_e_in)
            {
                size_t const &u = edge.first;
                size_t const &v = edge.second;
                ++k_node[u];
                ++k_node[v];
                ++M_in[u];
                ++M_in[v];

                number_of_edges++;
            }

            for(auto const &edge : *it_e_out)
            {
                size_t const &u = edge.first;
                size_t const &v = edge.second;
                --k_node[u];
                --k_node[v];
                ++M_out[u];
                ++M_out[v];

                number_of_edges--;
            }

            // advance network 
            ++it_t_ec;
            ++it_e_out;
            ++it_e_in;
        }

        if (t_next_ab <= t_next_ec)
        {
            // advance gamma and P
            ++it_g;
            ++it_P;
            ++it_t_ab;
        }

    } while ((t_next_ab < tmax) or (t_next_ec < tmax));


    vector < double > g_node;
    vector < double > P_node;
    vector < double > a_node;
    vector < double > b_node;


    for (size_t node = 0; node < N; ++node)
    {
        if (use_event_rate_method)
        {
            double this_a = (M_in[node] - I_in_1[node]) / I_in_2[node];
            double this_b = (M_out[node] - I_out[node] - I_out_2[node] - this_a * I_out[node]) / I_out_2[node];

            if (this_a < 0.0)
                this_a = 0.0;

            if (this_b < 0.0)
                this_b = 0.0;

            a_node.push_back(this_a);
            b_node.push_back(this_b);
        } else
        {
            double this_g = (M_out[node] - I_out[node]) / I_out[node];
            double this_P = (M_in[node] - I_in_1[node]) / I_in_2[node] / this_g;

            g_node.push_back(this_g);
            P_node.push_back(this_P);
        }
    }


    if (use_event_rate_method)
        return make_pair(a_node, b_node);
    else
        return make_pair(g_node, P_node);

}

flockwork_alpha_beta_args
     get_flockwork_alpha_beta_args(
             edge_changes &list_of_edge_changes,
             double dt,
             size_t N_time_steps,
             double k_over_k_real_scaling,
             double alpha_scaling,
             double beta_scaling,
             map < pair < size_t, size_t >, double > &aggregated_network,
             const bool ensure_empty_network,
             const bool adjust_last_bin_if_dt_does_not_fit,
             const bool use_integral_method,
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
    double original_tmax = tmax;

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
            if (adjust_last_bin_if_dt_does_not_fit)
            {
                // if it did not, change the time of the network
                // to be larger than originally.
                // This will be changed back later when 
                // computing gamma and P
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
    vector < double > M;

    size_t number_of_edges = list_of_edge_changes.edges_initial.size();


    new_time.push_back( t0 );

    if (tmax < time.back())
        throw domain_error("The value tmax is smaller than the last value in the time list.");

    if( verbose)
    {
        cout << "starting binning process with dt = " << dt << " and N_time_steps = " << N_time_steps << endl;
    }

    double last_time = t0;

    // the demanded new time bin is next one new_time = t0 + dt
    for(size_t t_i = 1; t_i < N_time_steps+1; t_i++)
    {
        double this_new_time = dt * t_i + t0;

        // initialize the counts for this time bin
        m_in.push_back(0);
        m_out.push_back(0);
        M.push_back( 0. );

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

            if (verbose)
            {
                cout << "edges going in = " << it_edges_in->size() << endl;
                cout << "edges going out = " << it_edges_out->size() << endl;
            }

            m_in.back() += it_edges_in->size();
            m_out.back() += it_edges_out->size();

            M.back() += (*it_time - last_time) * number_of_edges;
            
            number_of_edges += (it_edges_in->size() - it_edges_out->size());

            if (verbose)
            {
                cout << "number of edges = " << number_of_edges << endl;
                cout << "dt = " << (*it_time - last_time) << endl;
                cout << "M.back() = " << M.back() << endl;
            }

            last_time = *it_time;
            it_time++;
            it_edges_in++;
            it_edges_out++;
        }

        if ( (this_new_time == tmax) and original_tmax < tmax)
            this_new_time = original_tmax;

        M.back() += (this_new_time - last_time) * number_of_edges;

        new_time.push_back( this_new_time );
    }


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

    m_in.back() += incoming_edge_integers.size();
    m_out.back() += outgoing_edge_integers.size();

    //
    if (verbose)
    {
        cout << "the number of demanded time steps was " << N_time_steps << endl;
        cout << "the number of covered times in `new_time` was " << new_time.size() << endl << "    new_time = [ ";
        for(auto const & this_time: new_time)
        {            
            cout << this_time << ", ";
        }
        cout << "]" << endl;
        cout << "the number of covered M values in `M` was " << M.size() << endl << "    M = [ ";
        for(auto const & _M: M)
        {            
            cout << _M << ", ";
        }
        cout << "]" << endl;
        cout << "the number of covered m_in values in `m_in` was " << m_in.size() << endl << "    m_in = [ ";
        for(auto const & _m_in: m_in)
        {            
            cout << _m_in << ", ";
        }
        cout << "]" << endl;
        cout << "the number of covered m_out values in `m_out` was " << m_out.size() << endl << "    m_out = [ ";
        for(auto const & _m_out: m_out)
        {            
            cout << _m_out << ", ";
        }
        cout << "]" << endl;
    }

    // ========= now go through the observables and compute the parameters =============
    auto it_M = M.begin();
    auto it_m_in = m_in.begin();
    auto it_m_out = m_out.begin();
    it_time = new_time.begin();

    //initialize parameters
    vector < pair < double, double > > alpha;
    vector < double > beta;

    bool last_bin_was_for_emptying = false;

    while(it_m_in != m_in.end())
    {
        double & this_M = *it_M;
        size_t & this_m_in = *it_m_in;
        size_t & this_m_out = *it_m_out;

        // if tmax was shifted because dt did not fit, rescale the last
        // time dt to the appropriate dt
        if (it_m_in + 1 == m_in.end())
        {
            double this_time = *it_time;
            double next_time = original_tmax;
            dt = next_time - this_time;
        }


        if ( (this_M == 0.0) and ( this_m_out > 0) )
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

        double _a;
        double _b;

        double _m_out = (double) this_m_out;
        double _m_in = (double) this_m_in;
        double _M = this_M / k_over_k_real_scaling;

        _a = _m_in / (N*dt+2.0*_M);

        if ( (_m_out == 0.0) and (_M == 0.0) )
            _b = 0.0;
        else
            _b = (_m_out / 2.0 / _M) - _a;

        if (_b < 0.0)
            _b = 0.0;

        // if in the beginning of this bin and at the end of this bin
        // there's no edges and we want this to be ensured
        if ( ( ensure_empty_network ) and 
             ( this_M == 0.0 ) and 
             ( *(it_M+1) == 0.0 ) and
             ( _a == 0.0 ) and 
             ( _b == 0.0 )
           )
        {
            // we want to pick every node in the network at least once, s.t.
            // all nodes have zero edges in the end of the time bin.
            // This is the classic problem of "how often do I have to pick
            // in order to have picked every number at least once" with solution
            // <T> = N log(N)
            // Since we want <Events([t,t+dt])> = N * log(N), we need
            //  gamma_per_node = (N log(N)) / (N / dt) such that
            //  \int_t^{t+dt} dt' gamma_per_node * N = N * log(N)
            //
            //  Also, instead of doing this for two time bins, we double the rate
            //  and just do it for one time bin
            if ( not last_bin_was_for_emptying )
                _b = 1 * log(N) / dt;

            last_bin_was_for_emptying = true;
        }
        else
        {
            last_bin_was_for_emptying = false;
        }

        alpha.push_back( make_pair( *it_time, _a ) );
        beta.push_back( _b );

        it_M++; 
        it_m_in++;
        it_m_out++; 
        it_time++;
    }

    flockwork_alpha_beta_args fw_args;

    fw_args.N = N;
    fw_args.E = list_of_edge_changes.edges_initial;
    fw_args.disconnection_rate = beta;
    fw_args.reconnection_rate = alpha;
    fw_args.tmax = list_of_edge_changes.tmax;
    fw_args.new_time = new_time;
    fw_args.m_in = m_in;
    fw_args.m_out = m_out;
    fw_args.M = M;

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

pair < vector < double >, vector < double > >
    get_node_alpha_and_beta(
                edge_changes &ec,
                vector < pair < double, double > > &alpha,
                vector < double > &beta,
                double k_over_k_real_scaling,
                bool apply_mean_correction,
                bool verbose
                )
{
    // get references to edge_list and time
    vector < vector < pair < size_t, size_t > > > & all_edges_in = ec.edges_in;
    vector < vector < pair < size_t, size_t > > > & all_edges_out = ec.edges_out;
    vector < double > &ec_time = ec.t;

    // get vectors from gamma and P
    vector < double > ab_time;
    vector < double > a;

    for(auto &a_entry: alpha)
    {
        ab_time.push_back(a_entry.first);
        a.push_back(a_entry.second);
    }

    // set initial and final time
    double tmax = ec.tmax;

    // graph observables
    size_t N = ec.N;
    double Nf(N);
    size_t number_of_edges = ec.edges_initial.size();

    // intialize the node properties
    vector < size_t > k_node(N);
    vector < size_t > M_out(N);
    vector < size_t > M_in(N);

    // intialize the integrals
    vector < double > I_out_1(N);
    vector < double > I_out_2(N);
    vector < double > I_in_1(N);
    vector < double > I_in_2(N);

    for(auto const &edge : ec.edges_initial)
    {
        ++k_node[edge.first];
        ++k_node[edge.second];
    }

    auto it_e_in = all_edges_in.begin();
    auto it_e_out = all_edges_out.begin();

    auto it_t_ab = ab_time.begin();
    auto it_t_ec = ec_time.begin();

    auto it_a = a.begin();
    auto it_b = beta.begin();

    double t = ec.t0;

    double t_next_ab;
    double t_next_ec;
    double t_next;

    if (verbose)
        cout << "number of nodes = " << N << endl;

    do
    {

        // check which is the next event time at which stuff changes
        if ((it_t_ab == ab_time.end()) or (it_t_ab+1 == ab_time.end()))
            t_next_ab = tmax;
        else
            t_next_ab = *(it_t_ab+1);

        if ((it_t_ec == ec_time.end()) or (it_t_ec+1 == ec_time.end()))
            t_next_ec = tmax;
        else
            t_next_ec = *(it_t_ec+1);

        if ( t_next_ec <= t_next_ab )
            t_next = t_next_ec;
        else
            t_next = t_next_ab;

        // compute the integrals up to the next change

        double dt = t_next - t;
        double mean_degree = ((double) number_of_edges) * 2.0 / Nf / k_over_k_real_scaling;

        auto i_k = k_node.begin();
        auto i_I_out_1 = I_out_1.begin();
        auto i_I_out_2 = I_out_2.begin();
        auto i_I_in_1 = I_in_1.begin();
        auto i_I_in_2 = I_in_2.begin();

        double &_a = *it_a;
        double &_b = *it_b;

        //cout << "_a = " << _a << endl;
        //cout << "_b = " << _b << endl;

        while(i_k != k_node.end())
        {
            double _k = ((double) *i_k) / k_over_k_real_scaling;

            (*i_I_out_1) += _a * _k * dt;
            (*i_I_out_2) += _b * _k * dt;
            (*i_I_in_1) += _a * (Nf - _k - 1.0) * (_k + 1.0) / (Nf - 1.0) * dt;
            (*i_I_in_2) += _a * (1.0 + mean_degree) * dt;

            ++i_k;
            ++i_I_out_1;
            ++i_I_out_2;
            ++i_I_in_1;
            ++i_I_in_2;
        }

        // advance time
        t = t_next;

        // change either the network or gamma and P or both

        if (t_next_ec <= t_next_ab)
        {
            // update observables for this advancement
            for(auto const &edge : *it_e_in)
            {
                size_t const &u = edge.first;
                size_t const &v = edge.second;
                ++k_node[u];
                ++k_node[v];
                ++M_in[u];
                ++M_in[v];

                number_of_edges++;
            }

            for(auto const &edge : *it_e_out)
            {
                size_t const &u = edge.first;
                size_t const &v = edge.second;
                --k_node[u];
                --k_node[v];
                ++M_out[u];
                ++M_out[v];

                number_of_edges--;
            }

            // advance network 
            ++it_t_ec;
            ++it_e_out;
            ++it_e_in;
        }

        if (t_next_ab <= t_next_ec)
        {
            // advance <alpha> and <beta>
            ++it_a;
            ++it_b;
            ++it_t_ab;
        }

    } while ((t_next_ab < tmax) or (t_next_ec < tmax));


    vector < double > a_node;
    vector < double > b_node;

    double Alpha = 0.0;
    double Beta = 0.0;

    double mean_M_in = 0.0;
    double mean_M_out = 0.0;
    double mean_I_in_1 = 0.0;
    double mean_I_in_2 = 0.0;
    double mean_I_out_1 = 0.0;
    double mean_I_out_2 = 0.0;

    for (size_t node = 0; node < N; ++node)
    {
        if (verbose)
        {
            cout << "node = " << node << "; " << "M_in = " << M_in[node] << endl;
            cout << "node = " << node << "; " << "M_out = " << M_out[node] << endl;
            cout << "node = " << node << "; " << "I_in_1 = " << I_in_1[node] << endl;
            cout << "node = " << node << "; " << "I_out_1 = " << I_out_1[node] << endl;
            cout << "node = " << node << "; " << "I_in_2 = " << I_in_2[node] << endl;
            cout << "node = " << node << "; " << "I_out_2 = " << I_out_2[node] << endl;
        }
        mean_M_in += M_in[node] / (double) N;
        mean_M_out += M_out[node] / (double) N;
        mean_I_in_1 += I_in_1[node] / (double) N;
        mean_I_out_1 += I_out_1[node] / (double) N;
        mean_I_in_2 += I_in_2[node] / (double) N;
        mean_I_out_2 += I_out_2[node] / (double) N;
    }

    //double mean_alpha = mean_M_in / ( mean_I_in_1 + mean_I_in_2 );
    //double mean_beta = mean_M_out / ( mean_I_out_1 + mean_I_out_2 ) - mean_alpha;
    double mean_alpha = mean_M_in / ( mean_I_in_1 + mean_I_in_2 );
    double mean_beta = (mean_M_out*0.5 - mean_alpha) / mean_I_out_2;


    for (size_t node = 0; node < N; ++node)
    {

        /*
        cout << "M_in[" << node << "] = " << M_in[node] << endl;
        cout << "M_out[" << node << "] = " << M_out[node] << endl;
        cout << "I_in_1[" << node << "] = " << I_in_1[node] << endl;
        cout << "I_in_2[" << node << "] = " << I_in_2[node] << endl;
        cout << "I_out[" << node << "] = " << I_out[node] << endl;
        cout << "I_out_2[" << node << "] = " << I_out_2[node] << endl;
        */

        double this_a, this_b;

        if (not apply_mean_correction)
        {
            if (I_in_2[node] == 0.0)
                this_a = 0.0;
            else
                this_a = (M_in[node] - I_in_1[node]) / I_in_2[node];

            if (I_out_2[node] == 0.0)
                this_b = 0.0;
            else
                this_b = (M_out[node] - I_out_1[node] - I_out_2[node] - this_a * I_out_1[node]) / I_out_2[node];
        }
        else
        {
            if (I_in_2[node] == 0.0)
                this_a = 0.0;
            else
                this_a = (M_in[node] - mean_alpha * I_in_1[node]) / I_in_2[node];

            if (I_out_2[node] == 0.0)
                this_b = 0.0;
            else
                this_b = (M_out[node] - mean_alpha * I_out_1[node] - mean_beta * I_out_2[node] - this_a * I_out_1[node]) / I_out_2[node];
        }

        if (verbose)
        {
            cout << "node " << node << " ; alpha = "<< this_a << endl;
            cout << "node " << node << " ; beta  = "<< this_b << endl;
        }

        if (this_a < 0.0)
            this_a = 0.0;

        if (this_b < 0.0)
            this_b = 0.0;

        a_node.push_back(this_a);
        b_node.push_back(this_b);

        Alpha += this_a;
        Beta += this_b;

        if (verbose)
        {
            cout << "node " << node << " ; corrected alpha = "<< this_a << endl;
            cout << "node " << node << " ; corrected beta  = "<< this_b << endl;
        }

    }

    Alpha /= (double) N;
    Beta /= (double) N;

    for (size_t node = 0; node < N; ++node)
    {
        if (Alpha > 0.0)
            a_node[node] /= Alpha;
        if (Beta > 0.0)
            b_node[node] /= Beta;
    }

    return make_pair(a_node, b_node);
}

pair < vector < vector < double > >, vector < vector < double > > >
    get_time_dependent_node_alpha_and_beta(
                edge_changes &ec,
                vector < pair < double, double > > &alpha,
                vector < double > &beta,
                bool apply_mean_correction
                )
{
    // get references to edge_list and time
    vector < vector < pair < size_t, size_t > > > & all_edges_in = ec.edges_in;
    vector < vector < pair < size_t, size_t > > > & all_edges_out = ec.edges_out;
    vector < double > &ec_time = ec.t;

    // get vectors from gamma and P
    vector < double > ab_time;
    vector < double > a;

    for(auto &a_entry: alpha)
    {
        ab_time.push_back(a_entry.first);
        a.push_back(a_entry.second);
    }

    // set initial and final time
    double tmax = ec.tmax;

    // graph observables
    size_t N = ec.N;
    double Nf(N);
    size_t number_of_edges = ec.edges_initial.size();

    // intialize the node properties
    vector < size_t > k_node(N);
    vector < size_t > M_out(N);
    vector < size_t > M_in(N);

    // intialize the integrals
    vector < double > I_out_1(N);
    vector < double > I_out_2(N);
    vector < double > I_in_1(N);
    vector < double > I_in_2(N);

    for(auto const &edge : ec.edges_initial)
    {
        ++k_node[edge.first];
        ++k_node[edge.second];
    }

    auto it_e_in = all_edges_in.begin();
    auto it_e_out = all_edges_out.begin();

    auto it_t_ab = ab_time.begin();
    auto it_t_ec = ec_time.begin();

    auto it_a = a.begin();
    auto it_b = beta.begin();

    double t = ec.t0;

    double t_next_ab;
    double t_next_ec;
    double t_next;

    vector < vector < double > > alpha_t;
    vector < vector < double > > beta_t;

    do
    {

        // check which is the next event time at which stuff changes
        if ((it_t_ab == ab_time.end()) or (it_t_ab+1 == ab_time.end()))
            t_next_ab = tmax;
        else
            t_next_ab = *(it_t_ab+1);

        if ((it_t_ec == ec_time.end()) or (it_t_ec+1 == ec_time.end()))
            t_next_ec = tmax;
        else
            t_next_ec = *(it_t_ec+1);

        if ( t_next_ec <= t_next_ab )
            t_next = t_next_ec;
        else
            t_next = t_next_ab;

        // compute the integrals up to the next change

        double dt = t_next - t;
        double mean_degree = ((double) number_of_edges) * 2.0 / Nf;

        auto i_k = k_node.begin();
        auto i_I_out_1 = I_out_1.begin();
        auto i_I_out_2 = I_out_2.begin();
        auto i_I_in_1 = I_in_1.begin();
        auto i_I_in_2 = I_in_2.begin();

        double &_a = *it_a;
        double &_b = *it_b;

        //cout << "_a = " << _a << endl;
        //cout << "_b = " << _b << endl;

        while(i_k != k_node.end())
        {
            double _k = (double) *i_k;

            (*i_I_out_1) += _a * _k * dt;
            (*i_I_out_2) += _b * _k * dt;
            (*i_I_in_1) += _a * (Nf - _k - 1.0) * (_k + 1.0) / (Nf - 1.0) * dt;
            (*i_I_in_2) += _a * (1.0 + mean_degree) * dt;

            ++i_k;
            ++i_I_out_1;
            ++i_I_out_2;
            ++i_I_in_1;
            ++i_I_in_2;
        }

        // advance time
        t = t_next;

        // change either the network or gamma and P or both

        if (t_next_ec <= t_next_ab)
        {
            // update observables for this advancement
            for(auto const &edge : *it_e_in)
            {
                size_t const &u = edge.first;
                size_t const &v = edge.second;
                ++k_node[u];
                ++k_node[v];
                ++M_in[u];
                ++M_in[v];

                number_of_edges++;
            }

            for(auto const &edge : *it_e_out)
            {
                size_t const &u = edge.first;
                size_t const &v = edge.second;
                --k_node[u];
                --k_node[v];
                ++M_out[u];
                ++M_out[v];

                number_of_edges--;
            }

            // advance network 
            ++it_t_ec;
            ++it_e_out;
            ++it_e_in;
        }

        if (t_next_ab <= t_next_ec)
        {
            vector < double > a_node;
            vector < double > b_node;

            double Alpha = 0.0;
            double Beta = 0.0;

            double mean_M_in = 0.0;
            double mean_M_out = 0.0;
            double mean_I_in_1 = 0.0;
            double mean_I_in_2 = 0.0;
            double mean_I_out_1 = 0.0;
            double mean_I_out_2 = 0.0;

            for (size_t node = 0; node < N; ++node)
            {
                mean_M_in += M_in[node] / (double) N;
                mean_M_out += M_out[node] / (double) N;
                mean_I_in_1 += I_in_1[node] / (double) N;
                mean_I_out_1 += I_out_1[node] / (double) N;
                mean_I_in_2 += I_in_2[node] / (double) N;
                mean_I_out_2 += I_out_2[node] / (double) N;
            }

            //double mean_alpha = mean_M_in / ( mean_I_in_1 + mean_I_in_2 );
            //double mean_beta = mean_M_out / ( mean_I_out_1 + mean_I_out_2 ) - mean_alpha;
            double mean_alpha = mean_M_in / ( mean_I_in_1 + mean_I_in_2 );
            double mean_beta = (mean_M_out*0.5 - mean_alpha) / mean_I_out_2;


            if (( mean_I_in_1 + mean_I_in_2 ) == 0.0)
                mean_alpha = 0.0;

            if (( mean_I_out_1 + mean_I_out_2 ) == 0.0)
                mean_beta = 0.0;

            for (size_t node = 0; node < N; ++node)
            {

                double this_a, this_b;

                if (not apply_mean_correction)
                {
                    if (I_in_2[node] == 0.0)
                        this_a = 0.0;
                    else
                        this_a = (M_in[node] - I_in_1[node]) / I_in_2[node];

                    if (I_out_2[node] == 0.0)
                        this_b = 0.0;
                    else
                        this_b = (M_out[node] - I_out_1[node] - I_out_2[node] - this_a * I_out_1[node]) / I_out_2[node];
                }
                else
                {
                    if (I_in_2[node] == 0.0)
                        this_a = 0.0;
                    else
                        this_a = (M_in[node] - mean_alpha * I_in_1[node]) / I_in_2[node];

                    if (I_out_2[node] == 0.0)
                        this_b = 0.0;
                    else
                        this_b = (M_out[node] - mean_alpha * I_out_1[node] - mean_beta * I_out_2[node] - this_a * I_out_1[node]) / I_out_2[node];
                }

                if (this_a < 0.0)
                    this_a = 0.0;

                if (this_b < 0.0)
                    this_b = 0.0;

                a_node.push_back(this_a);
                b_node.push_back(this_b);

                Alpha += this_a;
                Beta += this_b;
            }

            Alpha /= (double) N;
            Beta /= (double) N;

            for (size_t node = 0; node < N; ++node)
            {
                if (Alpha > 0.0)
                    a_node[node] *= (*it_a) / Alpha;
                if (Beta > 0.0)
                    b_node[node] *= (*it_b) / Beta;

                // reset the integrals
                M_out[node] = 0.0;
                M_in[node] = 0.0;

                I_out_1[node] = 0.0;
                I_out_2[node] = 0.0;
                I_in_1[node] = 0.0;
                I_in_2[node] = 0.0;
            }

            alpha_t.push_back(a_node);
            beta_t.push_back(b_node);
            
            // advance gamma and P
            ++it_a;
            ++it_b;
            ++it_t_ab;

        }

    } while ((t_next_ab < tmax) or (t_next_ec < tmax));

    return make_pair(alpha_t, beta_t);
}
