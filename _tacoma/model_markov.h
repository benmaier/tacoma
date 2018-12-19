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

#ifndef __MODEL_MARKOV_H__
#define __MODEL_MARKOV_H__

#include "Utilities.h"

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
#include <assert.h>

using namespace std;

template <typename T1, typename T2>
void 
    markov_on_model(
            T1 & this_model_object,
            T2 & this_markov_object,
            double max_dt,
            bool reset_simulation_objects = true,
            bool verbose = false
            )
{

    if (this_model_object.N != this_markov_object.N)
        throw domain_error("Both model and markov object need to have the same node number.");

    if (verbose)
        cout << "started markov integration on model." << endl;

    // deal with random numbers
    mt19937_64 &generator = this_markov_object.generator;
    uniform_real_distribution<double> randuni(0.0,1.0);
    this_model_object.set_generator(this_markov_object.generator);


    // reset the simulation objects
    if (reset_simulation_objects)
    {
        if (verbose)
            cout << "Resetting models..." << endl;

        this_markov_object.reset();
        this_model_object.reset();

        if (verbose)
            cout << "Done." << endl;
    }

    // initialize time variables
    double t0 = this_model_object.t0;
    double t = t0;

    double t_simulation = this_markov_object.t_simulation;
    this_model_object.edg_chg.tmax = t_simulation;

    // initialize a graph and pass it to the markov object
    this_markov_object.update_network(this_model_object.G,t);

    if (verbose)
        cout << "Start integration" << endl;

    while ( (t-t0 < t_simulation) and (not this_markov_object.simulation_ended()) )
    {

        // inititalize rate containers
        vector < double > rates_model;
        double Lambda_model;

        // get the updated rates and lambda
        this_model_object.get_rates_and_Lambda(rates_model,Lambda_model);

        if (verbose)
        {
            cout << "Lambda_model = " << Lambda_model << endl;
            for(auto const &rate: rates_model)
                cout << "   rate = " << rate << endl;
        }

        double rProduct = randuni(generator) * Lambda_model;

        vector<double>::iterator this_rate;
        size_t n_rates;

        this_rate = rates_model.begin();
        n_rates = rates_model.size();

        double sum_event = 0.0;
        size_t event = 0;

        while ( (event < n_rates) and not ( (sum_event < rProduct) and (rProduct <= sum_event + (*this_rate)) ) )
        {
            sum_event += (*this_rate);
            ++this_rate;
            ++event;
        }

        if (verbose)
        {
            cout << "rProduct = " << rProduct << endl;
            cout << "event = " << event << endl;
            cout << "this_rate = " << *this_rate << endl;
        }

        double tau = log(1.0/randuni(generator)) / Lambda_model;

        while ((tau > max_dt) and (t-t0+max_dt < t_simulation) and (not this_markov_object.simulation_ended()) )
        {
            this_markov_object.step(t,max_dt);
            tau -= max_dt;
            t += max_dt;
        }

        if ((t-t0+tau < t_simulation) and (not this_markov_object.simulation_ended()))
        {
            this_markov_object.step(t, tau);

            vector < pair < size_t, size_t > > edges_in;
            vector < pair < size_t, size_t > > edges_out;

            this_model_object.make_event(event,t,edges_in,edges_out);
            this_markov_object.update_network(this_model_object.G,t+tau);
        }
        t += tau;

        if (verbose)
            cout << "new time: " << t << " / " << t_simulation << endl;

    }
}

#endif
