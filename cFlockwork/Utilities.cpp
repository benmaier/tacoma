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

#include "Utilities.h"


using namespace std;

pair <size_t,size_t> get_sorted_pair(size_t i, size_t j)
{
    //swap i and j if j is the smaller number
    if (j<i) 
        swap(i,j);

    return make_pair(i,j);
}

void choose (const size_t N, size_t &first, size_t &second, const double r1, const double r2)
{

  first = N * r1;
  second = (N-1) * r2;

  if (second >= first)
     second++;
}

//This function is based on http://ideone.com/3A3cv
//or http://stackoverflow.com/questions/9345087/choose-m-elements-randomly-from-a-vector-containing-n-elements

/*
template<class bidiiter> 
bidiiter choose_random_unique(
                 bidiiter begin, 
                 bidiiter end, 
                 size_t num_random, 
                 default_random_engine & generator, 
                 uniform_real_distribution<double> & distribution
               )
{
    size_t left = distance(begin, end);
    while (num_random--) {
        bidiiter r = begin;
        advance(r, distribution(generator)*left);
        swap(*begin, *r);
        ++begin;
        --left;
    }
    return begin;
}
*/

vector<size_t>::iterator choose_random_unique(
        vector<size_t>::iterator begin, 
        vector<size_t>::iterator end, 
        size_t num_random,
        default_random_engine & generator,
        uniform_real_distribution<double> & distribution
    ) 
{
    size_t left = distance(begin, end);
    while (num_random--) {
        vector<size_t>::iterator r = begin;
        advance(r, distribution(generator) * left);
        swap(*begin, *r);
        ++begin;
        --left;
    }
    return begin;
}

size_t arg_choose_from_vector(
        vector < double > const & weights, 
        default_random_engine & generator, 
        uniform_real_distribution<double> & distribution
      )
{
    double a0 = accumulate(weights.begin(), weights.end(), 0.0);
    double rProduct = distribution(generator) * a0;

    int event = 0;
    int N = weights.size();

    if (N==0)
    {
        throw length_error( "Rate list is empty." );
    }

    double sum_event = 0.0;
    while ( (event<N) and not ( (sum_event < rProduct) and (rProduct <= sum_event+weights[event]) ) )
    {
        sum_event += weights[event];
        event++;
    }

    return event;
}


void get_gillespie_tau_and_event( 
                vector < double > const & rates,
                double & tau,
                size_t & event,
                default_random_engine & generator, 
                uniform_real_distribution<double> & distribution
             )
{
    double a0 = accumulate(rates.begin(), rates.end(), 0.0);
    double rProduct = distribution(generator) * a0;

    event = 0;
    size_t N = rates.size();

    if (N==0)
    {
        throw length_error( "Rate list is empty." );
    }

    double sum_event = 0.0;
    while ( (event<N) and not ( (sum_event < rProduct) and (rProduct <= sum_event+rates[event]) ) )
    {
        sum_event += rates[event];
        event++;
    }

    tau = log(1.0/distribution(generator)) / a0;
}

void get_gillespie_tau_and_event_with_varying_gamma( 
                vector < double > & standard_rates,
                vector < pair < double, double > > const & gamma,
                double t0,
                double t_max,
                size_t & i_t,
                double & tau,
                size_t & event,
                default_random_engine & generator, 
                uniform_real_distribution<double> & distribution
             )
{
    // ======================= FIND TAU ========================
    size_t N_gamma = gamma.size();
    double beta0 = accumulate(standard_rates.begin(), standard_rates.end(), 0.0);
    double m_log_1mU = - log( 1. - distribution(generator));

    double t_basic = (i_t / N_gamma) * t_max;

    double t_iP1 = get<0>(gamma[(i_t+1) % N_gamma]) + ((i_t+1) / N_gamma) * t_max; 
    double t_i   = t0;
    double g_i   = get<1>(gamma[i_t % N_gamma]);

    double S_i = 0.;
    double S_iP1 = (t_iP1 - t0) * g_i;

    double next_limit = (t_iP1 - t0) * beta0 + S_iP1;

    while (m_log_1mU > next_limit) 
    {
        i_t++;
        t_basic = (i_t / N_gamma) * t_max;

        g_i   = get<1>(gamma[i_t % N_gamma]);
        t_i   = get<0>(gamma[i_t % N_gamma]) + t_basic;
        t_iP1 = get<0>(gamma[(i_t+1) % N_gamma]) + ((i_t+1) / N_gamma) * t_max;

        S_i = S_iP1;
        S_iP1 += (t_iP1 - t_i) * g_i;
        next_limit = (t_iP1 - t0) * beta0 + S_iP1;
    }


    //tau = ( m_log_1mU - S_i + g_i*t_i + t0*beta0 ) / (beta0 + g_i) - t0;
    //
    tau = (m_log_1mU - S_i + g_i * (t_i-t0) ) / (beta0 + g_i);
    double Lambda = g_i * (tau+(t0-t_i)) + S_i;

    if ((t0+tau) == t_iP1)
        i_t++;


    //================ FIND EVENT ========================
    size_t N = standard_rates.size();
    double a0 = Lambda;

    //get the mean number of events happened per channel
    for(size_t rate=0; rate<N; rate++)
    {
        standard_rates[rate] *= tau;
        a0 += standard_rates[rate];
    }

    //push rewiring event to the end
    standard_rates.push_back(Lambda);
    N++;

    //find event as usual
    double rProduct = distribution(generator) * a0;

    event = 0;

    if (N==0)
    {
        throw length_error( "Rate list is empty." );
    }

    double sum_event = 0.0;
    while ( (event<N) and not ( (sum_event < rProduct) and (rProduct <= sum_event+standard_rates[event]) ) )
    {
        sum_event += standard_rates[event];
        event++;
    }

}

