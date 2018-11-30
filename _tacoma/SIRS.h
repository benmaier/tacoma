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

#ifndef __SIRS_H__
#define __SIRS_H__

#include "Utilities.h"
#include "ResultClasses.h"

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

class SIRS 
{
    public:
        size_t N;
        double t_simulation;
        double infection_rate;
        double recovery_rate;
        double becoming_susceptible_rate;
        size_t number_of_initially_infected;
        size_t number_of_initially_vaccinated;
        size_t seed;
        bool verbose;
        double sampling_dt;

        mt19937_64 generator;
        uniform_real_distribution<double> randuni;

        vector < double > time;
        vector < double > R0;
        vector < size_t > SI;
        vector < size_t > I;
        vector < size_t > R;

        SIRS(
            size_t _N,
            double _t_simulation,
            double _infection_rate,
            double _recovery_rate,
            double _becoming_susceptible_rate,
            size_t _number_of_initially_infected = 1,
            size_t _number_of_initially_vaccinated = 0,
            double _sampling_dt = 0.0,
            size_t _seed = 0,
            bool _verbose = false
        )
        {
            N = _N;
            t_simulation = _t_simulation;
            infection_rate = _infection_rate;
            recovery_rate = _recovery_rate;
            becoming_susceptible_rate = _becoming_susceptible_rate;
            number_of_initially_vaccinated = _number_of_initially_vaccinated;
            number_of_initially_infected = _number_of_initially_infected;
            verbose = _verbose;
            seed = _seed;
            sampling_dt = _sampling_dt;

            mt19937_64 generator;

            randuni = uniform_real_distribution < double > (0.0, 1.0);
        }

        void reset() 
        {
            next_sampling_time = sampling_dt;

            // reset observables
            //
            time.clear();
            R0.clear();
            SI.clear();
            I.clear();
            R.clear();

            // seed engine
            if (seed == 0)
                randomly_seed_engine(generator);
            else
                generator.seed(seed);

            //check if number of infected and number of vaccinated does not
            //exceed total node number
            if (number_of_initially_vaccinated + number_of_initially_infected > N) 
                throw length_error( "Number of infected and number of vaccinated may not exceed total population size" );

            //initialize status vector of nodes and vector of infected
            node_status = vector < size_t >(N,EPI::S);
            infected = vector < size_t >();
            recovered = vector < size_t >();

            vector < size_t > node_ints;
            for(size_t n=0; n<N; n++)
            {
                node_ints.push_back(n);            
            }

            if (verbose) 
            {
                cout << "choosing " << number_of_initially_vaccinated
                     << " vaccinated and " << number_of_initially_infected
                     << " infected at random" << endl;
            }

            choose_random_unique(node_ints.begin(),
                                 node_ints.end(),
                                 number_of_initially_vaccinated + number_of_initially_infected,
                                 generator,
                                 randuni
                                 );

            for(size_t node=0; node<number_of_initially_vaccinated; node++)
            {
                node_status[node_ints[node]] = EPI::V;
            }

            for(size_t node = number_of_initially_vaccinated; node < number_of_initially_vaccinated + number_of_initially_infected; node++)
            {
                node_status[node_ints[node]] = EPI::I;
                infected.push_back( node_ints[node] );
            }

            if (verbose)
            {
                cout << "infected set has size = " << infected.size() << endl;
                print_infected();
            }
        }

        bool simulation_ended() 
        {
            return ( (recovered.size() == 0) and (infected.size() == 0) );
        };

        void get_rates_and_Lambda(vector < double > &rates,
                                  double &Lambda
                                  );

        void update_network(vector < set < size_t > > &G,
                            double t
                            );

        void update_network(vector < set < size_t > > &G,
                            vector < pair < size_t, size_t > > &edges_in,
                            vector < pair < size_t, size_t > > &edges_out,
                            double t
                            );

        void make_event(size_t const &event,
                        double t
                       );

        void update_observables(double t);

        void print()
        {
            print_infected();
            print_SI_edges();
        }

    private:
        vector < size_t > infected;
        vector < size_t > recovered;
        vector < size_t > node_status;
        vector < pair < size_t, size_t > > SI_edges;
        double mean_degree;
        vector < set < size_t > > * G;

        double next_sampling_time;

        vector < double > rates;

        void infection_event();
        void recovery_event();
        void susceptible_event();

        void print_recovered()
        {
            cout << "recovered = [ ";
            for( auto const &rec: recovered)
                cout << rec << " ";
            cout << "]" << endl;
        };

        void print_infected()
        {
            cout << "infected = [ ";
            for( auto const &inf: infected)
                cout << inf << " ";
            cout << "]" << endl;
        };

        void print_SI_edges()
        {
            cout << "SI_edges = [ ";
            for( auto const &edge: SI_edges)
                cout << "( "<< edge.first << ", " << edge.second << " ) ";
            cout << "]" << endl;
        };

};

#endif
