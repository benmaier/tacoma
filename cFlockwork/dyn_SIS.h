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

#ifndef __DYN_SIS_H__
#define __DYN_SIS_H__

#include "Utilities.h"
#include "Events.h"

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

const size_t _S = 0;
const size_t _I = 1;
const size_t _R = 1;

struct Dyn_SIS 
{
    size_t N;
    double t_simulation;
    double infection_rate;
    double recovery_rate;
    size_t number_of_initially_infected;
    size_t number_of_initially_vaccinated;
    size_t seed;

    default_random_engine generator;
    uniform_real_distribution<double> randuni;

    vector < double > t;
    vector < double > R0;
    vector < size_t > SI;
    vector < size_t > I;

    set < size_t > infected;
    vector < size_t > node_status;
    set < pair < size_t, size_t > > SI_E;

    vector < double > rates;

    Dyn_SIS(
        size_t _N,
        double _t_simulation,
        double _infection_rate,
        double _recovery_rate,
        size_t _number_of_initially_infected = 1,
        size_t _number_of_initially_vaccinated = 0,
        size_t _seed = 0
    )
    {
        N = _N;
        t_simulation = _t_simulation;
        infection_rate = _infection_rate;
        recovery_rate = _recovery_rate;
        number_of_initially_vaccinated = _number_of_initially_vaccinated;
        number_of_initially_infected = _number_of_initially_infected;

        default_random_engine generator;

        if (_seed == 0)
            randomly_seed_engine(generator);
        else
            generator.seed(_seed);

        randuni = uniform_real_distribution<double>(0.0,1.0);

        //check if number of infected and number of vaccinated does not
        //exceed total node number
        if (number_of_initially_vaccinated + number_of_initially_infected > N) 
            throw length_error( "Number of infected and number of vaccinated may not exceed total population size" );

        //initialize status vector of nodes and vector of infected
        vector < size_t > node_status(N,_S);
        set < size_t > infected;


        vector < size_t > node_ints;
        for(size_t n=0; n<N; n++)
        {
            node_ints.push_back(n);            
        }

        choose_random_unique(node_ints.begin(),
                             node_ints.end(),
                             number_of_initially_vaccinated + number_of_initially_infected,
                             generator,
                             randuni
                             );

        for(size_t node=0; node<number_of_initially_vaccinated; node++)
        {
            node_status[node_ints[node]] = _R;
        }

        for(size_t node = number_of_initially_vaccinated; node < number_of_initially_vaccinated + number_of_initially_infected; node++)
        {
            node_status[node_ints[node]] = _I;
            infected.insert( node_ints[node] );
        }


    }


    bool simulation_ended() {
        return infected.size() == 0;
    };

    void update_network(vector < set < size_t > > G,
                        double t
                        ){};
    void get_rates_and_Lambda(vector < double > &rates,
                              double &Lambda
                              ){};
    void make_event(size_t const &event,
                    double t
                   ){};

    void update_observables(){};
};


#endif
