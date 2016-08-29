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
#include "SIS.h"

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


        
//returns vector containing vaccinated 
tuple < 
        vector < pair < double, size_t > >, 
        vector < pair < double, size_t > >, 
        vector < pair < double, double > >,
        set < pair < size_t, size_t >  >
      >
     SIS(
                 vector < tuple < size_t, size_t > > E, //edgelist
                 const size_t N,       //number of nodes
                 const double Q,       //probability to connect with neighbors of neighbor
                 const size_t t_run_total,
                 const double infection_rate,
                 const double recovery_rate,
                 const double rewiring_rate,
                 const size_t number_of_vaccinated,
                 const size_t number_of_infected,
                 const size_t seed
        )
{

    //check if number of infected and number of vaccinated does not
    //exceed total node number
    if (number_of_vaccinated + number_of_infected > N) 
        throw length_error( "Number of infected and number of vaccinated may not exceed total population size" );

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

    for(size_t node=0; node<N; node++)
        G.push_back(new set < size_t >);

    //initialize edge list of infected-susceptible links
    set < pair < size_t, size_t > > SI_E;

    //loop through edge list and push neighbors
    for(auto edge: E)
    {
        //get nodes belonging to that edge
        size_t i = get<0>(edge);
        size_t j = get<1>(edge);

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

    //initialize random generators
    default_random_engine generator(seed);
    uniform_real_distribution<double> uni_distribution(0.,1.);

    //init result vectors 
    vector < pair <double,size_t> > I_of_t;
    vector < pair <double,size_t> > SI_of_t;
    vector < pair <double,double> > R0_of_t;

    I_of_t.push_back(make_pair(0.,infected.size()));
    SI_of_t.push_back(make_pair(0.,SI_E.size()));
    R0_of_t.push_back(make_pair(0.,k*infection_rate/recovery_rate));

    
    //simulate
    double t = 0;
    size_t last_event = -1;
    while ( (t < t_run_total) && (infected.size()>0) )
    {
        //calculate rates
        vector <double> rates;
        rates.push_back(N*rewiring_rate);
        rates.push_back(SI_E.size()*infection_rate);
        rates.push_back(infected.size()*recovery_rate);

        double tau;
        size_t event;
        get_gillespie_tau_and_event(rates,tau,event,generator,uni_distribution);
        t = t + tau;
        last_event = event;

        if (event == 0)
        {
            rewire(G,Q,generator,uni_distribution,k,SI_E,node_status);
            R0_of_t.push_back(make_pair(t,k*infection_rate/recovery_rate));
        } else if (event == 1)
        {
            infect(G,generator,uni_distribution,SI_E,node_status,infected);
            I_of_t.push_back(make_pair(t,infected.size()));
            SI_of_t.push_back(make_pair(t,SI_E.size()));
        } else if (event == 2)
        {
            SIS_recover(G,generator,uni_distribution,SI_E,node_status,infected);
            I_of_t.push_back(make_pair(t,infected.size()));
            SI_of_t.push_back(make_pair(t,SI_E.size()));
        }
    }

    if (last_event == 0)
    {
        I_of_t.push_back( make_pair( t, infected.size() ) );
        SI_of_t.push_back( make_pair( t, SI_E.size() ) );
    } else {
        R0_of_t.push_back(make_pair(t,k*infection_rate/recovery_rate));
    }


    //convert back to edge list
    set < pair < size_t, size_t > > new_E;    

    //loop through graph
    for(size_t node=0; node<N; node++)
    {

        //loop through nodes' neighbors and add edge
        for(auto const& neigh: *G[node])
            new_E.insert( get_sorted_pair(node,neigh)  );

        //delete created set
        delete G[node];
    }

    return make_tuple(I_of_t, SI_of_t, R0_of_t, new_E);
}

