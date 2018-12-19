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

#include "FlockworkPModel.h"

using namespace std;
//namespace osis = SIS;
void
    FlockworkPModel::rewire(
                 double const &_P,       //probability to connect with neighbors of neighbor
                 vector < pair < size_t, size_t > > &edges_in,
                 vector < pair < size_t, size_t > > &edges_out
            )
{
    //choose two nodes
    double r1 = randuni(*generator);
    double r2 = randuni(*generator);
    size_t i,j;
    choose(N,i,j,r1,r2);

    bool do_rewiring = randuni(*generator) < _P;

    //check if new neighbor is actually an old neighbor
    //and if this is the case return an empty event
    if (not ( do_rewiring and (G[i].find(j) != G[i].end()) ))
    {
        //loop through the neighbors of i
        for(auto neigh_i : G[i] )
        {
            //and erase the link to i
            G[neigh_i].erase(i);
            
            pair < size_t, size_t > current_edge = get_sorted_pair(i,neigh_i);
            edges_out.push_back( current_edge );
        } 

        //erase the links from the perspective of i
        G[i].clear();

        if ( do_rewiring )
        {
            //loop through the neighbors of j
            for(auto neigh_j : G[j] ) 
            {
                G[ neigh_j ].insert( i );
                G[ i ].insert( neigh_j );

                pair < size_t, size_t > current_edge = get_sorted_pair(i,neigh_j);
                edges_in.push_back( current_edge );
            }

            //add j as neighbor and count additional edge
            G[ j ].insert( i );
            G[ i ].insert( j );

            pair < size_t, size_t > current_edge = get_sorted_pair(i,j);
            edges_in.push_back( current_edge );
        }
    }
}


void FlockworkPModel::get_rates_and_Lambda(
                    vector < double > &_rates,
                    double &_Lambda
                  )
{
    // delete current rates
    _rates = rates;
    _Lambda = Lambda;
}

void FlockworkPModel::make_event(
                size_t const &event,
                double t,
                vector < pair < size_t, size_t > > &e_in,
                vector < pair < size_t, size_t > > &e_out
               )
{

    e_in.clear();
    e_out.clear();

    if (event==0)
        rewire(1.0,e_in,e_out);
    else if (event == 1)
        rewire(0.0,e_in,e_out);
}
