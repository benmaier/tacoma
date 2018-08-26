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

#include "Events.h"
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

edge_changes
     activity_model(
                 const size_t N,       //number of nodes
                 const double rho,       //probability to reconnect after cutting
                 const double omega,
                 const double t_run_total,
                 const size_t seed,
                 const bool   verbose
        )
{
    assert((rho>0.0) and (rho <1.0));
    assert(omega > 0.0);

    double omega_minus = omega / rho;
    double omega_plus = omega / (1.0 - rho);

    if (verbose)
    {
        cout << "omega = " << omega << endl;
        cout << "rho = " << rho << endl;
        cout << "omega+ = " << omega_plus << endl;
        cout << "omega- = " << omega_minus << endl;
    }

    //initialize random generators
    mt19937_64 generator;
    seed_engine(generator,seed);
    uniform_real_distribution<double> uni_distribution(0.,1.);

    vector < size_t > k(N);
    vector < size_t > complementary_k(N);

    vector < set < size_t > > G = get_random_graph(N, rho, generator);

    for(size_t node=0; node<N; node++)
    {
        k[node] = G[node].size();
        complementary_k[node] = N - 1 - k[node];

        if (verbose)
        {
            cout << "k[node]   = " << k[node] << endl;
            cout << "c_k[node] = " << complementary_k[node] << endl;
        }
    }
    
    vector < pair < size_t, size_t > > edges_initial;
    edgelist_from_graph(edges_initial, G);

    vector < vector < pair <size_t,size_t> > > edges_out;
    vector < vector < pair <size_t,size_t> > > edges_in;
    vector < double > time;

    double t_initial = 0.0;
    double t = t_initial;

    size_t edges_on = edges_initial.size();
    size_t m_max = (N * (N-1)) / 2;

    while (t < t_run_total)
    {
        //calculate rates
        vector <double> rates;
        rates.push_back(edges_on * omega_minus);
        rates.push_back((m_max - edges_on) * omega_plus);

        double tau;
        size_t event;
        get_gillespie_tau_and_event(rates, tau, event, generator, uni_distribution);

        if (verbose)
        {
            cout << "=========================" << endl;
            for(auto const &rate: rates)
                cout << "rate = " << rate << endl;
            cout << "tau = " << tau << endl;
            cout << "event = " << event << endl;
            cout << endl;
        }

        t += tau;

        if (t<t_run_total)
        {
            vector < pair < size_t, size_t > > e_out;
            vector < pair < size_t, size_t > > e_in;

            if (event==0)
            {
                // find a turned on edge and turn it off.
                
                // find a reacting node (probability proportional to degree)
                size_t reacting_node = arg_choose_from_vector(k, generator, uni_distribution);

                if (verbose)
                {
                    cout << "found node " << reacting_node << " with degree " << k[reacting_node] << endl;
                    cout << "neighs = ";
                    for(auto const &neigh : G[reacting_node])
                        cout << neigh << " ";
                    cout << endl;
                }

                // find a random neighbor of this node and iterate to it
                uniform_int_distribution < size_t > random_neighbor(0,k[reacting_node]-1);

                auto it = G[reacting_node].begin();

                for(size_t i = 0; i < random_neighbor(generator); i++)
                    ++it;

                size_t reacting_neighbor = *it;

                if (verbose)
                    cout << "found neighbor " << reacting_neighbor << endl;

                // add this edge to edges going out
                e_out.push_back(get_sorted_pair(reacting_node, reacting_neighbor));

                // adjust degree and complementary degree
                edges_on--;

                k[reacting_node]--;
                k[reacting_neighbor]--;
                complementary_k[reacting_node]++;
                complementary_k[reacting_neighbor]++;

                G[reacting_node].erase(reacting_neighbor);
                G[reacting_neighbor].erase(reacting_node);
            }
            else if (event == 1)
            {
                //find a turned off edge and turn it on.
                
                // find a reacting node (probability proportional to degree)
                size_t reacting_node = arg_choose_from_vector(complementary_k, generator, uni_distribution);
                set < size_t > &neighs =  G[reacting_node];

                if (verbose)
                {
                    cout << "found node " << reacting_node << " with complementary degree " << complementary_k[reacting_node] << endl;
                    cout << "neighs = ";
                    for(auto const &neigh : neighs)
                        cout << neigh << " ";
                    cout << endl;
                    cout << "inserting temporary self-link to avoid picking it" << endl;
                }

                neighs.insert(reacting_node);


                // find a random non-neighbor index of this node
                uniform_int_distribution < size_t > random_neighbor(0,complementary_k[reacting_node]-1);
                size_t reacting_neighbor = random_neighbor(generator);
                if (verbose)
                    cout << "   current reacting non-neighbor candidate = " << reacting_neighbor << endl;

                // find an index of a random non-neighbor
                auto it_neighbor = neighs.begin();

                // iterate through neighbors
                while ( (it_neighbor != neighs.end()) and 
                        (*it_neighbor <= reacting_neighbor)
                      )
                {
                    if (verbose)
                    {
                        cout << "   current neighbor = " << *it_neighbor << endl;
                        cout << "   current reacting non-neighbor candidate = " << reacting_neighbor << endl;
                    }
                    // while the index of the non-neighbor is larger than the index of the current_neighbor,
                    // increase the index of the non-neighbor by one and look at the next neighbor
                    reacting_neighbor++;
                    ++it_neighbor;

                }

                neighs.erase(reacting_node);

                if (verbose)
                    cout << "found reacting non-neighbor = " << reacting_neighbor << endl;

                // add this edge to edges going out
                e_in.push_back(get_sorted_pair(reacting_node, reacting_neighbor));

                // adjust degree and complementary degree
                edges_on++;

                k[reacting_node]++;
                k[reacting_neighbor]++;
                complementary_k[reacting_node]--;
                complementary_k[reacting_neighbor]--;

                G[reacting_node].insert(reacting_neighbor);
                G[reacting_neighbor].insert(reacting_node);
            }
            else
            {
                throw length_error("There was an event chosen other than on or off, this should not happen.");
            }

            edges_out.push_back(e_out);
            edges_in.push_back(e_in);
            time.push_back(t);

            if(verbose)
            {
                cout << "graph is now" << endl;
                for(size_t node = 0; node < N; node++)
                {
                    cout << "G[" << node << "] = ";
                    for(auto const &neigh: G[node])
                        cout << neigh << " ";
                    cout << endl;
                }
                cout << "degrees/complementary are now" <<endl;
                for(size_t node = 0; node < N; node++)
                {
                    cout << "  k[" << node << "]   = " << k[node] << endl;
                    cout << "  c_k[" << node << "] = " << complementary_k[node] << endl;
                }
            }

        }
    }

    edge_changes result;

    result.t = time;
    result.tmax = t_run_total;
    result.t0 = t_initial;
    result.edges_initial = edges_initial;
    result.edges_out = edges_out;
    result.edges_in = edges_in;
    result.N = N;

    return result;
}
