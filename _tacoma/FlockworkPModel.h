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

#ifndef __FLOCKWORKP_MODEL_CLASS_H__
#define __FLOCKWORKP_MODEL_CLASS_H__

#include "Utilities.h"
#include "ResultClasses.h"
#include "activity_model.h"

//#include <iostream>
//#include <algorithm>
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

class FlockworkPModel 
{
    public:
        size_t N;
        double gamma;
        double P;
        double alpha;
        double beta;
        double t0;
        size_t seed;
        bool verbose;
        bool save_temporal_network;
        vector < pair < size_t, size_t > > E;
        vector < double > rates;
        double Lambda;


        vector < set < size_t > > G;

        mt19937_64 * generator;
        uniform_real_distribution<double> randuni;

        edge_changes edg_chg;

        FlockworkPModel(
            vector < pair < size_t, size_t > > _E,
            size_t _N,
            double _gamma,
            double _P,
            double _t0 = 0.0, 
            bool _save_temporal_network = false,
            size_t _seed = 0,
            bool _verbose = false
        )
        {
            E = _E;
            N = _N;
            P = _P;
            gamma = _gamma;
            t0 = _t0;
            verbose = _verbose;
            seed = _seed;
            save_temporal_network = _save_temporal_network;

            alpha = gamma*P;
            beta = gamma*(1-P);

            rates.push_back(N*alpha);
            rates.push_back(N*beta);
            Lambda = N*alpha + N*beta;


            generator = new mt19937_64;
            randuni = uniform_real_distribution < double > (0.0, 1.0);

            if (verbose)
            {
                cout << "gamma = " << gamma << endl;
                cout << "P = " << P << endl;
                cout << "alpha = " << alpha << endl;
                cout << "beta = " << beta << endl;
            }

        }

        void set_initial_configuration( double _t0, vector < set < size_t > > &_G)
        {
            t0 = _t0;

            // reset observables
            edg_chg.N = N;

            edg_chg.t.clear();
            edg_chg.edges_in.clear();
            edg_chg.edges_out.clear();

            // seed engine
            if (seed == 0)
                randomly_seed_engine(*generator);
            else
                generator->seed(seed);


            G.clear();
            for(size_t node=0; node<N; node++)
            {
                set < size_t > neigh = _G[node];
                G.push_back(neigh);
            }

            if (save_temporal_network)
                edgelist_from_graph(edg_chg.edges_initial, G);

        }

        void reset() 
        {
            // reset observables
            edg_chg.N = N;

            edg_chg.t.clear();
            edg_chg.edges_in.clear();
            edg_chg.edges_out.clear();

            // seed engine
            if (seed == 0)
                randomly_seed_engine(*generator);
            else
                generator->seed(seed);

            G = vector < set < size_t > >(N);
            graph_from_edgelist(G, E);


            if (save_temporal_network)
                edg_chg.edges_initial = E;
        }

        void get_rates_and_Lambda(vector < double > &rates,
                                  double &Lambda
                                  );

        void make_event(size_t const &event,
                        double t,
                        vector < pair < size_t, size_t > > &e_in,
                        vector < pair < size_t, size_t > > &e_out
                       );

        void print()
        {
        }

        void set_generator(mt19937_64 &_generator)
        {
            delete generator;
            generator = &_generator;
            has_external_generator = true;
        } 

        ~FlockworkPModel() 
        {
            if (not has_external_generator)
                delete generator;
        }

    private:

        bool has_external_generator = false;

        void
            rewire(
                     double const &_P,       //probability to connect with neighbors of neighbor
                     vector < pair < size_t, size_t > > &edges_in,
                     vector < pair < size_t, size_t > > &edges_out
                  );

        void print_internal()
        {
        }



};

#endif
