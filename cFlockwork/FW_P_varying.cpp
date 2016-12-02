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
#include "SIS_varying.h"
#include "EqFlockwork.h"

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

const size_t S = 0;
const size_t I = 1;
const size_t R = 2;
const size_t E_OUT = 0;
const size_t E_IN = 1;


        
edge_changes
     flockwork_P_varying_rates(
                 vector < pair < size_t, size_t > > E, //edgelist
                 const size_t N,       //number of nodes
                 vector < double > P,       //probability to reconnect after cutting
                 const double t_run_total,
                 vector < pair < double, double > > rewiring_rate,
                 const double tmax,
                 const bool   use_random_rewiring,
                 const bool   equilibrate_flockwork,
                 const size_t seed
        )
{

    double infection_rate = 0;
    double recovery_rate = 0;
    size_t number_of_vaccinated = 0;
    size_t number_of_infected = 1;

    //check if number of infected and number of vaccinated does not
    //exceed total node number
    if (number_of_vaccinated + number_of_infected > N) 
        throw length_error( "Number of infected and number of vaccinated may not exceed total population size" );

    //initialize random generators
    default_random_engine generator(seed);
    uniform_real_distribution<double> uni_distribution(0.,1.);

    if (equilibrate_flockwork)
    {
        size_t eq_time = 0;
        E = equilibrate_edgelist_generator(E,N,P[0],generator,uni_distribution,eq_time,true);
    }

    //initialize status vector of nodes and vector of infected
    vector < size_t > node_status;
    set < size_t > infected;

    for(size_t node=0; node<number_of_vaccinated; node++)
        node_status.push_back( R );

    for(size_t node=number_of_vaccinated; node<number_of_vaccinated+number_of_infected; node++)
    {
        node_status.push_back( I );
        infected.insert( node );
    }

    for(size_t node=number_of_vaccinated+number_of_infected; node<N; node++)
        node_status.push_back( S );

    //initialize Graph vector
    vector < set < size_t > * > G;

    //count number of edges in Graph
    size_t number_of_edges = 0;

    //if we use random rewiring, we need a vector containing the node ints
    //this seems to be very bad style but I don't have a better idea right now
    vector < size_t > node_ints;

    for(size_t node=0; node<N; node++)
    {
        G.push_back(new set < size_t >);
        node_ints.push_back(node);
    }

    //initialize edge list of infected-susceptible links
    set < pair < size_t, size_t > > SI_E;

    //loop through edge list and push neighbors
    for(auto edge: E)
    {
        //get nodes belonging to that edge
        size_t i = edge.first;
        size_t j = edge.second;

        //check if j is already a neighbor of i
        const bool already_counted = G[i]->find(j) != G[i]->end();

        //check if edge has been added already
        if (!(already_counted)) 
        {
            G[ i ]->insert( j );
            G[ j ]->insert( i );
            number_of_edges++;

            //check if this is an SI link
            if ( 
                 ( (node_status[i] == I) && (node_status[j] == S) ) ||
                 ( (node_status[i] == S) && (node_status[j] == I) )
               ) 
            {
                pair <size_t,size_t> current_pair = get_sorted_pair(i,j);
                SI_E.insert( current_pair );
            }
        }

    }

    //calculate mean degree
    double k = (2.0 / N) * number_of_edges;

    vector < vector < pair <size_t,size_t> > > edges_out;
    vector < vector < pair <size_t,size_t> > > edges_in;
    vector < double > time;

    //multiply rewiring rate with number of nodes
    size_t N_gamma = rewiring_rate.size();
    for (size_t it = 0; it<N_gamma; it++)
        get<1>(rewiring_rate[it]) *= N;

    if ( N_gamma != P.size() )
        throw length_error( "Q and rewiring rate need to have the same length." );

    //simulate
    double t = get<0>(rewiring_rate[0]);
    size_t last_event = -1;
    size_t i_t = 0;

    while ( (t < t_run_total) && (infected.size()>0) )
    {
        //calculate rates
        vector <double> rates;
        rates.push_back(infected.size()*recovery_rate);

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
        last_event = event;

        if (event==1)
        {
            pair < vector < pair <size_t,size_t> >, 
                   vector < pair <size_t,size_t> > 
                 > curr_edge_change = rewire_P(G,P[i_t%N_gamma],generator,uni_distribution,k,SI_E,node_status);

            if ( (curr_edge_change.first.size()>0) || (curr_edge_change.second.size()>0) )
            {
                time.push_back(t);
                edges_out.push_back( curr_edge_change.first );
                edges_in.push_back( curr_edge_change.second );
            }
            /*

            cout << t << endl;

            cout << "out" << endl;

            for(auto edge: curr_edge_change.first)
                cout << edge.first << " " << edge.second<< endl;

            cout << "in" << endl;
            for(auto edge: curr_edge_change.second)
                cout << edge.first << " " << edge.second<< endl;

            for(size_t u=0; u<N; u++)
                for( auto v: *G[u])
                {
                    cout << u << " " << v << endl;

                }
            */
        }
        else
        {
            throw length_error("There was an event chosen other than rewiring, this should not happen.");
        }
            

    }

    edge_changes result;

    result.t = time;
    result.edges_out = edges_out;
    result.edges_in = edges_in;

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
                 const bool   equilibrate_flockwork,
                 const bool   use_preferential_node_selection,
                 const bool   use_unweighted_k_for_selection,
                 const size_t seed
        )
{

    double infection_rate = 0;
    double recovery_rate = 0;
    size_t number_of_vaccinated = 0;
    size_t number_of_infected = 1;

    //check if number of infected and number of vaccinated does not
    //exceed total node number
    if (number_of_vaccinated + number_of_infected > N) 
        throw length_error( "Number of infected and number of vaccinated may not exceed total population size" );

    //initialize random generators
    default_random_engine generator(seed);
    uniform_real_distribution<double> uni_distribution(0.,1.);

    if (equilibrate_flockwork)
    {
        size_t eq_time = 0;
        E = equilibrate_edgelist_generator(E,N,P[0],generator,uni_distribution,eq_time,true);
    }

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

    //initialize status vector of nodes and vector of infected
    vector < size_t > node_status;
    set < size_t > infected;

    for(size_t node=0; node<number_of_vaccinated; node++)
        node_status.push_back( R );

    for(size_t node=number_of_vaccinated; node<number_of_vaccinated+number_of_infected; node++)
    {
        node_status.push_back( I );
        infected.insert( node );
    }

    for(size_t node=number_of_vaccinated+number_of_infected; node<N; node++)
        node_status.push_back( S );

    //initialize Graph vector
    vector < set < size_t > * > G;

    //count number of edges in Graph
    size_t number_of_edges = 0;

    //if we use random rewiring, we need a vector containing the node ints
    //this seems to be very bad style but I don't have a better idea right now
    vector < size_t > node_ints;

    for(size_t node=0; node<N; node++)
    {
        G.push_back(new set < size_t >);
        node_ints.push_back(node);
    }

    //initialize edge list of infected-susceptible links
    set < pair < size_t, size_t > > SI_E;

    //loop through edge list and push neighbors
    for(auto edge: E)
    {
        //get nodes belonging to that edge
        size_t i = edge.first;
        size_t j = edge.second;

        //check if j is already a neighbor of i
        const bool already_counted = G[i]->find(j) != G[i]->end();

        //check if edge has been added already
        if (!(already_counted)) 
        {
            G[ i ]->insert( j );
            G[ j ]->insert( i );
            number_of_edges++;

            //check if this is an SI link
            if ( 
                 ( (node_status[i] == I) && (node_status[j] == S) ) ||
                 ( (node_status[i] == S) && (node_status[j] == I) )
               ) 
            {
                pair <size_t,size_t> current_pair = get_sorted_pair(i,j);
                SI_E.insert( current_pair );
            }
        }

    }

    //calculate mean degree
    double k = (2.0 / N) * number_of_edges;

    vector < vector < pair <size_t,size_t> > > edges_out;
    vector < vector < pair <size_t,size_t> > > edges_in;
    vector < double > time;

    //multiply rewiring rate with number of nodes
    size_t N_gamma = rewiring_rate.size();
    for (size_t it = 0; it<N_gamma; it++)
        get<1>(rewiring_rate[it]) *= N;

    if ( N_gamma != P.size() )
        throw length_error( "Q and rewiring rate need to have the same length." );

    //simulate
    double t = get<0>(rewiring_rate[0]);
    size_t last_event = -1;
    size_t i_t = 0;

    while ( (t < t_run_total) && (infected.size()>0) )
    {
        //calculate rates
        vector <double> rates;
        rates.push_back(infected.size()*recovery_rate);

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
        last_event = event;

        if (event==1)
        {
            pair < vector < pair <size_t,size_t> >, 
                   vector < pair <size_t,size_t> > 
                 > curr_edge_change = rewire_P_neighbor_affinity(G,
                                                                 P[i_t%N_gamma],
                                                                 neighbor_affinity,
                                                                 total_affinity,
                                                                 generator,
                                                                 uni_distribution,
                                                                 k,
                                                                 SI_E,
                                                                 node_status,
                                                                 use_preferential_node_selection);

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

    edge_changes result;

    result.t = time;
    result.edges_out = edges_out;
    result.edges_in = edges_in;

    return result;
}

