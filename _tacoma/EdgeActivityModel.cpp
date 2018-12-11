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

#include "EdgeActivityModel.h"

using namespace std;
//namespace sis = SIS;

void EdgeActivityModel::get_rates_and_Lambda(
                    vector < double > &_rates,
                    double &_Lambda
                  )
{
    // delete current rates
    _rates.clear();

    if (verbose)
    {
        cout << "edges_on = " << edges_on << endl;
        cout << "omega- = " << omega_minus << endl;
        cout << "omega+ = " << omega_plus << endl;
        cout << "OMEGA- = " << omega_minus * ((double) edges_on) << endl; 
        cout << "OMEGA+ = " << omega_plus * ((double) (m_max - edges_on)) << endl;
    }

    // compute rates of infection
    _rates.push_back(omega_minus * ((double) edges_on));
    _rates.push_back(omega_plus * ((double) (m_max-edges_on)));

    // return those new rates
    _Lambda = accumulate(_rates.begin(),_rates.end(),0.0);
}

void EdgeActivityModel::make_event(
                size_t const &event,
                double t,
                vector < pair < size_t, size_t > > &e_in,
                vector < pair < size_t, size_t > > &e_out
               )
{

    e_in.clear();
    e_out.clear();

    if (event==0)
    {
        // find a turned on edge and turn it off.
        
        // find a reacting node (probability proportional to degree)
        size_t reacting_node = arg_choose_from_vector(k, *generator, randuni);

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
        advance(it, random_neighbor(*generator));

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
        size_t reacting_node = arg_choose_from_vector(complementary_k, *generator, randuni);
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
        size_t reacting_neighbor = random_neighbor(*generator);
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
        throw length_error("EdgeActivityModel: There was an event chosen other than on or off, this should not happen.");
    }

    if (save_temporal_network)
    {
        edg_chg.edges_out.push_back(e_out);
        edg_chg.edges_in.push_back(e_in);
        edg_chg.t.push_back(t);
    }

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
