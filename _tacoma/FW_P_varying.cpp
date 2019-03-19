/* 
 * The MIT License (MIT)
 * Copyright (c) 2016, Benjamin Maier
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

#include "Events.h"
#include "Utilities.h"
#include "ResultClasses.h"
#include "Flockwork.h"
#include "FW_P_varying.h"

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
#include <ctime>
#include <tuple>

using namespace std;

//const size_t S = 0;
//const size_t I = 1;
//const size_t R = 2;

        
edge_changes
     flockwork_P_varying_rates(
                 vector < pair < size_t, size_t > > &E, //edgelist
                 const size_t N,       //number of nodes
                 vector < double > &P,       //probability to reconnect after cutting
                 const double t_run_total,
                 vector < pair < double, double > > &rewiring_rate,
                 const double tmax,
                 const bool   use_random_rewiring,
                 const size_t seed
        )
{

    //initialize random generators
    mt19937_64 generator;
    seed_engine(generator,seed);
    uniform_real_distribution<double> uni_distribution(0.,1.);

    //initialize Graph vector
    vector < set < size_t > * > G;

    //if we use random rewiring, we need a vector containing the node ints
    //this seems to be very bad style but I don't have a better idea right now
    vector < size_t > node_ints;

    for(size_t node=0; node<N; node++)
    {
        G.push_back(new set < size_t >);
        node_ints.push_back(node);
    }

    graph_from_edgelist(G, E);

    vector < vector < pair <size_t,size_t> > > edges_out;
    vector < vector < pair <size_t,size_t> > > edges_in;
    vector < double > time;

    //multiply rewiring rate with number of nodes
    size_t N_gamma = rewiring_rate.size();
    for (size_t it = 0; it<N_gamma; it++)
        get<1>(rewiring_rate[it]) *= N;

    if ( N_gamma != P.size() )
        throw length_error( "P and rewiring rate need to have the same length." );

    //simulate
    double t = get<0>(rewiring_rate[0]);
    size_t i_t = 0;

    while (t < t_run_total)
    {
        //calculate rates
        vector <double> rates;
        rates.push_back(0.0);

        double tau;
        size_t event;
        get_gillespie_tau_and_event_with_varying_gamma(
                            rates,
                            rewiring_rate,
                            t,
                            tmax,
                            i_t,
                            tau,
                            event,
                            generator,
                            uni_distribution);
        t = t + tau;

        if (t<t_run_total)
        {

            if (event==1)
            {
                pair < vector < pair <size_t,size_t> >, 
                       vector < pair <size_t,size_t> > 
                     > curr_edge_change = rewire_P_without_SI_checking(G,P[i_t%N_gamma],generator,uni_distribution);

                if ( (curr_edge_change.first.size()>0) || (curr_edge_change.second.size()>0) )
                {
                    time.push_back(t);
                    edges_out.push_back( curr_edge_change.first );
                    edges_in.push_back( curr_edge_change.second );
                }
            }
            else
            {
                throw length_error("There was an event chosen other than rewiring, this should not happen.");
            }
        }
            

    }

    edge_changes result;

    result.t = time;
    result.tmax = t_run_total;
    result.t0 = 0;
    result.edges_initial = E;
    result.edges_out = edges_out;
    result.edges_in = edges_in;
    result.N = N;

    return result;
}

edge_changes
     flockwork_P_varying_rates_for_each_node(
                 vector < pair < size_t, size_t > > &E, //edgelist
                 const size_t N,       //number of nodes
                 vector < vector < double > > &Ps,       //probability to reconnect after cutting
                 const double t_run_total,
                 vector < pair < double, vector < double > > > &rewiring_rates,
                 const double tmax,
                 const bool   use_random_rewiring,
                 const size_t seed
        )
{

    //initialize random generators
    mt19937_64 generator;
    seed_engine(generator,seed);
    uniform_real_distribution<double> uni_distribution(0.,1.);

    //initialize Graph vector
    vector < set < size_t > * > G;

    //if we use random rewiring, we need a vector containing the node ints
    //this seems to be very bad style but I don't have a better idea right now
    vector < size_t > node_ints;

    for(size_t node=0; node<N; node++)
    {
        G.push_back(new set < size_t >);
        node_ints.push_back(node);
    }

    graph_from_edgelist(G, E);

    vector < vector < pair <size_t,size_t> > > edges_out;
    vector < vector < pair <size_t,size_t> > > edges_in;
    vector < double > time;

    //multiply rewiring rate with number of nodes
    size_t N_gamma = rewiring_rates.size();
    vector < pair < double, double > > rewiring_rate_sum;

    for(auto &gamma_entry : rewiring_rates)
    {
        double this_gamma = accumulate(gamma_entry.second.begin(), gamma_entry.second.end(), 0.0);
        rewiring_rate_sum.push_back(make_pair(gamma_entry.first, this_gamma));
    }

    if ( N_gamma != Ps.size() )
        throw length_error( "P and rewiring rate need to have the same length." );

    //simulate
    double t = get<0>(rewiring_rates[0]);
    size_t i_t = 0;

    while (t < t_run_total)
    {
        //calculate rates
        vector <double> rates;
        rates.push_back(0.0);

        double tau;
        size_t event;

        //cout << "finding event using Gillespie SSA... " << endl;

        get_gillespie_tau_and_event_with_varying_gamma_for_each_node(
                            rates,
                            rewiring_rate_sum,
                            rewiring_rates,
                            t,
                            tmax,
                            i_t,
                            tau,
                            event,
                            generator,
                            uni_distribution);
        t = t + tau;

        if (t<t_run_total)
        {

            if (event>0)
            {
                size_t node = event - 1;

                //cout << "found node " << node << endl;
                //cout << "Ps[i_t%N_gamma].size() = " << Ps[i_t%N_gamma].size() << endl;
                
                pair < vector < pair <size_t,size_t> >, 
                       vector < pair <size_t,size_t> > 
                     > curr_edge_change = rewire_P_without_SI_checking_single_node(node,G,Ps[i_t%N_gamma][node],generator,uni_distribution);

                //cout << "rewired." << endl;

                if ( (curr_edge_change.first.size()>0) || (curr_edge_change.second.size()>0) )
                {
                    time.push_back(t);
                    edges_out.push_back( curr_edge_change.first );
                    edges_in.push_back( curr_edge_change.second );
                }
            }
            else
            {
                throw length_error("There was an event chosen other than rewiring, this should not happen.");
            }
        }
            

    }

    edge_changes result;

    result.t = time;
    result.tmax = t_run_total;
    result.t0 = 0;
    result.edges_initial = E;
    result.edges_out = edges_out;
    result.edges_in = edges_in;
    result.N = N;

    return result;
}

        
edge_changes
     flockwork_P_varying_rates_neighbor_affinity(
                 vector < pair < size_t, size_t > > E, //edgelist
                 const size_t N,       //number of nodes
                 vector < double > P,       //probability to reconnect after cutting
                 const double t_run_total,
                 vector < pair < double, double > > rewiring_rate,
                 vector < pair < vector < size_t >, vector < double > > > neighbor_affinity,
                 const double tmax,
                 const bool   use_random_rewiring,
                 const bool   use_preferential_node_selection,
                 const bool   use_unweighted_k_for_selection,
                 const size_t seed
        )
{

    //initialize random generators
    mt19937_64 generator;
    seed_engine(generator,seed);
    uniform_real_distribution<double> uni_distribution(0.,1.);

    vector < double > total_affinity(N);

    for(size_t node=0; node<N; node++)
    {
        double k_total;

        if (use_unweighted_k_for_selection)
            k_total = neighbor_affinity[node].first.size();
        else
            k_total = accumulate(neighbor_affinity[node].second.begin(), neighbor_affinity[node].second.end(), 0.0);

        total_affinity[node] = k_total;

    }

    //initialize Graph vector
    vector < set < size_t > * > G;

    //count number of edges in Graph
    //size_t number_of_edges = 0;

    //if we use random rewiring, we need a vector containing the node ints
    //this seems to be very bad style but I don't have a better idea right now
    vector < size_t > node_ints;

    for(size_t node=0; node<N; node++)
    {
        G.push_back(new set < size_t >);
        node_ints.push_back(node);
    }

    graph_from_edgelist(G, E);

    vector < vector < pair <size_t,size_t> > > edges_out;
    vector < vector < pair <size_t,size_t> > > edges_in;
    vector < double > time;

    //multiply rewiring rate with number of nodes
    size_t N_gamma = rewiring_rate.size();
    for (size_t it = 0; it<N_gamma; it++)
        get<1>(rewiring_rate[it]) *= N;

    if ( N_gamma != P.size() )
        throw length_error( "P and rewiring rate need to have the same length." );

    //simulate
    double t = get<0>(rewiring_rate[0]);
    size_t i_t = 0;

    while (t < t_run_total)
    {
        //calculate rates
        vector <double> rates;
        rates.push_back(0.0);

        double tau;
        size_t event;
        get_gillespie_tau_and_event_with_varying_gamma(
                            rates,
                            rewiring_rate,
                            t,
                            tmax,
                            i_t,
                            tau,
                            event,
                            generator,
                            uni_distribution);
        t = t + tau;

        if (t<t_run_total)
        {

            if (event==1)
            {
                pair < vector < pair <size_t,size_t> >, 
                       vector < pair <size_t,size_t> > 
                     > curr_edge_change = rewire_P_neighbor_affinity_without_SI_checking(
                                                                     G,
                                                                     P[i_t%N_gamma],
                                                                     neighbor_affinity,
                                                                     total_affinity,
                                                                     generator,
                                                                     uni_distribution,
                                                                     use_preferential_node_selection
                                                                     );

                if ( (curr_edge_change.first.size()>0) || (curr_edge_change.second.size()>0) )
                {
                    time.push_back(t);
                    edges_out.push_back( curr_edge_change.first );
                    edges_in.push_back( curr_edge_change.second );
                }
            }
            else
            {
                throw length_error("There was an event chosen other than rewiring, this should not happen.");
            }
        }
            

    }

    edge_changes result;

    result.t = time;
    result.tmax = t_run_total;
    result.t0 = 0;
    result.edges_initial = E;
    result.N = N;
    result.edges_out = edges_out;
    result.edges_in = edges_in;

    return result;
}

