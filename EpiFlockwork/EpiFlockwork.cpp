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
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace std;
namespace py = pybind11;

const size_t S = 0;
const size_t I = 1;
const size_t R = 2;
        

// ================================================ header =======================================

void choose (const size_t N, size_t &first, size_t &second, const double r1, const double r2);

pair <size_t,size_t> get_sorted_pair(size_t i, size_t j);

void rewire(
                 vector < set < size_t > * > & G, //Adjacency matrix
                 double Q,       //probability to connect with neighbors of neighbor
                 default_random_engine & generator, 
                 uniform_real_distribution<double> & distribution,
                 double & mean_degree,
                 set < pair < size_t, size_t > > & SI_E, //edge list of SI links
                 const vector < size_t > & node_status
              );

void infect(
                 vector < set < size_t > * > & G, //Adjacency matrix
                 default_random_engine & generator, 
                 uniform_real_distribution<double> & distribution,
                 set < pair < size_t, size_t > > & SI_E, //edge list of SI links
                 vector < size_t > & node_status,
                 set < size_t > & infected
           );

void SIS_recover(
                 vector < set < size_t > * > & G, //Adjacency matrix
                 default_random_engine & generator, 
                 uniform_real_distribution<double> & distribution,
                 set < pair < size_t, size_t > > & SI_E, //edge list of SI links
                 vector < size_t > & node_status,
                 set < size_t > & infected
           );

tuple < 
        vector < pair < double, size_t > >, 
        vector < pair < double, size_t > >, 
        vector < pair < double, double > >,
        vector < vector < size_t > * >
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
        );

// =============================================
//
void infect(
                 vector < set < size_t > * > & G, //Adjacency matrix
                 default_random_engine & generator, 
                 uniform_real_distribution<double> & distribution,
                 set < pair < size_t, size_t > > & SI_E, //edge list of SI links
                 vector < size_t > & node_status,
                 set < size_t > & infected
           )
{
    //get element of SI_E that leads to infection
    double r = distribution(generator);
    size_t element = SI_E.size() * r;

    set< pair < size_t, size_t > >::const_iterator edge(SI_E.begin());
    advance(edge,element);

    //get nodes belonging to that edge
    size_t i = (*edge).first;
    size_t j = (*edge).second;

    size_t new_infected;

    if ( ( node_status[i] == I ) && ( node_status[j] == S) )
    {
        new_infected = j;
    } else if ( ( node_status[j] == I ) && ( node_status[i] == S) )
    {
        new_infected = i;
    } else {
        throw domain_error( "There was a non SI link in the set of SI links. This should not happen." );
    }

    //change status of that node
    node_status[new_infected] = I;
    infected.insert(new_infected);

    //loop through the neighbors of i
    for(auto neigh_of_new_I : *G[new_infected] )
    {
        //get pair
        pair current_edge = get_sorted_pair(new_infected,neigh_of_new_I);

        //if the neighbor is infected, then the edge does not belong in SI_E anymore
        if (node_status[neigh_of_new_I] == I)
            SI_E.erase(current_edge);
        //if the neighbor is susceptible, then it belongs in SI_E now
        else if (node_status[neigh_of_new_I] == S)
            SI_E.insert(current_edge);
    }

}
                
void SIS_recover(
                 vector < set < size_t > * > & G, //Adjacency matrix
                 default_random_engine & generator, 
                 uniform_real_distribution<double> & distribution,
                 set < pair < size_t, size_t > > & SI_E, //edge list of SI links
                 vector < size_t > & node_status,
                 set < size_t > & infected
           )
{
    //get element of infected that recovers
    double r = distribution(generator);
    size_t element = infected.size() * r;

    set< size_t> >::const_iterator new_recovered(infected.begin());
    advance(new_recovered,element);

    //change status of that node
    node_status[new_recovered] = S;
    infected.erase(new_recovered);

    //loop through the neighbors of i
    for(auto neigh_of_new_S : *G[new_infected] )
    {
        //get pair
        pair current_edge = get_sorted_pair(new_infected,neigh_of_new_S);

        //if the neighbor is infected, then the edge belongs in SI_E
        if (node_status[neigh_of_new_S] == I)
            SI_E.insert(current_edge);
        //if the neighbor is susceptible, then it belongs in SI_E now
        else if (node_status[neigh_of_new_S] == S)
            SI_E.erase(current_edge);
    }

}

void rewire(
                 vector < set < size_t > * > & G, //Adjacency matrix
                 double Q,       //probability to connect with neighbors of neighbor
                 default_random_engine & generator, 
                 uniform_real_distribution<double> & distribution,
                 double & mean_degree,
                 set < pair < size_t, size_t > > & SI_E, //edge list of SI links
                 const vector < size_t > & node_status
            )
{
    //choose two nodes
    size_t N = G.size();
    double r1 = distribution(generator);
    double r2 = distribution(generator);
    size_t i,j;
    choose(N,i,j,r1,r2);

    size_t number_of_old_edges = 0;
    size_t umber_of_new_edges = 0;

    size_t number_of_old_SI_edges = 0;
    size_t umber_of_new_SI_edges = 0;

    //loop through the neighbors of i
    for(auto neigh_i : *G[i] )
    {
        //and erase the link to i
        G[neigh_i]->erase(i);
        number_of_old_edges++;

        //check if we erased an SI link
        if ( 
             ( (node_status[i] == I) && (node_status[neigh_i] == S) ) ||
             ( (node_status[i] == S) && (node_status[neigh_i] == I) )
           ) 
        {
            pair current_pair = get_sorted_pair(i,neigh_i);
            SI_E.erase( current_pair );
        }
    } 

    //erase the links from the perspective of i
    G[i]->clear();

    //loop through the neighbors of j
    for(auto neigh_j : *G[j] ) 
    {
        //add as neighbor of i if a random number is smaller than Q
        if ( distribution(generator) < Q )
        {
            G[ neigh_j ]->insert( i );
            G[ i ]->insert( neigh_j );
            number_of_new_edges++;

            //check if we created an SI link
            if ( 
                 ( (node_status[i] == I) && (node_status[neigh_j] == S) ) ||
                 ( (node_status[i] == S) && (node_status[neigh_j] == I) )
               ) 
            {
                pair current_pair = get_sorted_pair(i,neigh_j);
                SI_E.insert( current_pair );
            }
        }
    }

    //add j as neighbor and count additional edge
    G[ j ]->insert( i );
    G[ i ]->insert( j );
    number_of_new_edges++;

    //check if we created an SI link
    if ( 
         ( (node_status[i] == I) && (node_status[j] == S) ) ||
         ( (node_status[i] == S) && (node_status[j] == I) )
       ) 
    {
        pair current_pair = get_sorted_pair(i,j);
        SI_E.insert( current_pair );
    }

    //calculate new number of edges
    mean_degree += 2.0 * ( double(number_of_new_edges) - double(number_of_old_edges) ) / double(N);

}

pair <size_t,size_t> get_sorted_pair(size_t i, size_t j)
{
    //swap i and j if j is the smaller number
    if (j<i) 
    {
        size_t mem = i;
        i = j;
        j = mem;
    }

    return make_pair(i,j);
}


//returns vector containing vaccinated 
tuple < 
        vector < pair < double, size_t > >, 
        vector < pair < double, size_t > >, 
        vector < pair < double, double > >,
        vector < vector < size_t > * >
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

    t_max = 1000;
    //equilibrate
    for(size_t t=0; t<t_max; t++)
    {
        rewire(G,Q,generator,uni_distribution,k,SI_E,node_status);
    }

    //convert back to edge list
    vector < vector < size_t > * > new_E;


    //loop through graph
    for(size_t node=0; node<N; node++)
    {
        //loop through nodes' neighbors and add edge
        for(auto const& neigh: *G[node])
        {
            new_E.push_back(new vector < size_t >(2));
            new_E.back()->at(0) = node;
            new_E.back()->at(1) = neigh;
        }

        //delete created set
        delete G[node];
    }

    return new_E;
}


void choose (const size_t N, size_t &first, size_t &second, const double r1, const double r2)
{

  first = N * r1;
  second = (N-1) * r2;

  if (second >= first)
     second++;
}

PYBIND11_PLUGIN(EpiFlockwork) {
    py::module m("EpiFlockwork", "Module to equilibrate a flockwork in a fast manner");

    m.def("SIS", &SIS, "Simulate an SIS process on a flockwork given an initial state as an edge list. Returns time and number of infected as well as time and current R0.",
            py::arg("E"),
            py::arg("N"),
            py::arg("Q"),
            py::arg("t_run_total"),
            py::arg("infection_rate"),
            py::arg("recovery_rate"),
            py::arg("rewiring_rate") = 1.,
            py::arg("number_of_vaccinated") = 0,
            py::arg("number_of_infected") = 1,
            py::arg("seed") = 0,
            );

    return m.ptr();

}

int main()
{
    size_t N = 500;
    vector < tuple < size_t, size_t  > > E;
    double Q = 0.1;
    size_t t_run_total = 10;
    double infection_rate = 1.;
    double recovery_rate = 0.01;
    double rewiring_rate = 1.;
    size_t number_of_vaccinated = 0;
    size_t number_of_infected = 1;
    size_t seed = 1;

    tuple < 
            vector < pair < double, size_t > >, 
            vector < pair < double, size_t > >, 
            vector < pair < double, double > >,
            vector < vector < size_t > * >
          > result;

    result = SIS(E,N,Q,
                 t_run_total,
                 infection_rate,
                 recovery_rate,
                 rewiring_rate,
                 number_of_vaccinated,
                 number_of_infected,
                 seed
            );


    //for(size_t i=0; i<result.size(); i++)
    //{
    //    cout << result[i]->at(0) << " " << result[i]->at(1) << "\n";
    //}


    return 0;
}
