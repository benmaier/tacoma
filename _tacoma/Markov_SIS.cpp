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

#include "Markov_SIS.h"

using namespace std;
//namespace MARKOV_SIS = MARKOV_SIS;

void MARKOV_SIS::update_network(
                    vector < set < size_t > > &_G,
                    double t
                    )
{
    // map this network to the current network
    G = &_G;

    // compute degree and R0
    size_t number_of_edges_times_two = 0;

    for(auto const &neighbors: *G)
        number_of_edges_times_two += neighbors.size();

    mean_degree = number_of_edges_times_two / (double) N;

    // update the arrays containing the observables
    update_observables(t);
}

void MARKOV_SIS::step(
        double t,
        double dt
        )
{
    //perform step
    auto this_p_inf = p_infected->begin();
    auto status = node_status.begin();
    vector < double > *new_p = new vector < double >;

    double P_INF = dt * infection_rate;
    double P_REC = dt * recovery_rate;

    //cout << "P_INF = " << P_INF << endl;
    //cout << "P_REC = " << P_REC << endl;


    current_I = 0.0;
    for(auto const &neighbors: *G)
    {
        size_t const &stat = *status;
        double this_new_p = 0.0;
        if (stat != EPI::R)
        {
            double const &_p = *this_p_inf;
            double this_noninfection_probability = 1.0;
            for(auto const &j: neighbors)
            {
                if (node_status[j] != EPI::R)
                {
                    this_noninfection_probability *= 1.0 - P_INF * p_infected->at(j);
                }
            }
            this_new_p =   (1.0 - _p) * (1.0 - this_noninfection_probability)
                         + _p * (1 - P_REC);

            //cout << "node_status : " << stat << endl;
        }

        new_p->push_back(this_new_p);

        current_I += this_new_p;

        ++this_p_inf;
        ++status;
    }

    //cout << "current_I = " << current_I << endl;

    delete p_infected;
    p_infected = new_p;

    //update time and observables
    update_observables(t+dt);
}

void MARKOV_SIS::update_observables(
                double t
               )
{
    if (sampling_dt > 0.0)
    {
        if ((t >= next_sampling_time) and (time.back()!=t))
        {
            double _R0 = infection_rate * mean_degree / recovery_rate;
            R0.push_back(_R0);

            // compute I
            I.push_back(current_I);

            // push back time
            time.push_back(t);

            // advance next sampling time
            do
            {
                next_sampling_time += sampling_dt;
            } while (next_sampling_time < t);
        }
    }
    else
    {
        double _R0 = infection_rate * mean_degree / recovery_rate;
        R0.push_back(_R0);

        // compute I
        I.push_back(current_I);

        // push back time
        time.push_back(t);
    }
}

