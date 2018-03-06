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

#include "dyn_SIS.h"

using namespace std;
//namespace sis = Dyn_SIS;

void Dyn_SIS::update_network(
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

    size_t old_N_SI = SI_edges.size();

    // update infected
    SI_edges.clear();

    // for each infected, check its neighbors.
    // if the neighbor is susceptible, push it to the
    // endangered susceptible vector
    for(auto const &inf: infected)
        for(auto const &neighbor_of_infected: (*G)[inf])
            if (node_status[neighbor_of_infected] == _S)
                SI_edges.push_back( make_pair( inf, neighbor_of_infected ) );

    // update the arrays containing the observables
    update_observables(t);
}

void Dyn_SIS::get_rates_and_Lambda(
                    vector < double > &_rates,
                    double &_Lambda
                  )
{
    // delete current rates
    rates.clear();

    // compute rates of infection
    rates.push_back(infection_rate * SI_edges.size());

    // compute rates of recovery
    rates.push_back(recovery_rate * infected.size());

    // return those new rates
    _rates = rates;
    _Lambda = accumulate(rates.begin(),rates.end(),0.0);
}

void Dyn_SIS::make_event(
                size_t const &event,
                double t
               )
{
    if (event == 0)
        infection_event();
    else if (event == 1)
        recovery_event();
    else
        throw length_error("Dyn_SIS: chose event larger than rate vector which should not happen.");

    update_observables(t);
}

void Dyn_SIS::infection_event()
{
    // initialize uniform integer random distribution
    uniform_int_distribution < size_t > random_susceptible(0,SI_edges.size()-1);

    // find the index of the susceptible which will become infected
    size_t this_susceptible_index = random_susceptible(generator);

    if (verbose)
    {
        cout << "====================== INFECTION EVENT =====================" << endl;
        cout << "found SI edge with index " << this_susceptible_index << endl;
    }

    // get the node number of this susceptible
    size_t this_susceptible = (SI_edges.begin() + this_susceptible_index)->second;

    // save this node as an infected
    infected.push_back(this_susceptible);

    if (verbose)
        cout << "changing node status of susceptible " << this_susceptible << " to infected" << endl;

    // change node status of this node
    node_status[this_susceptible] = _I;

    if (verbose)
        cout << "erasing all SI edges which " << this_susceptible << " was part of" << endl;

    // erase all edges in the SI set where this susceptible is part of
    SI_edges.erase( 
            remove_if( 
                    SI_edges.begin(), 
                    SI_edges.end(),
                [&this_susceptible](const pair < size_t, size_t > & edge) { 
                        return edge.second == this_susceptible;
                }),
            SI_edges.end() 
        );

    if (verbose)
    {
        cout << " SI edges after erasing " << endl;
        print_SI_edges();
        cout << "finding all edges which " << this_susceptible << " is part of and is connected to a susceptible " << endl;
    }

    // push the new SI edges
    size_t & this_infected = this_susceptible;
    for(auto const &neighbor: (*G)[this_infected])
    {
        if (verbose)
            cout << "new infected " << this_susceptible << " has neighbor " << neighbor << " with status " << node_status[neighbor] << endl;

        if (node_status[neighbor] == _S)
            SI_edges.push_back( make_pair(this_infected, neighbor) );
    }

    if (verbose)
    {
        cout << "  SI edges after pushing = " << endl;
        print_SI_edges();
    }
}

void Dyn_SIS::recovery_event()
{
    // initialize uniform integer random distribution
    uniform_int_distribution < size_t > random_infected(0,infected.size()-1);

    // find the index of the susceptible which will become infected
    size_t this_infected_index = random_infected(generator);
    auto it_infected = infected.begin() + this_infected_index;

    if (verbose)
    {
        cout << "====================== RECOVERY EVENT =====================" << endl;
        cout << "found infected with index " << this_infected_index << endl;
    }

    // get the node id of this infected about to be recovered
    size_t this_infected = *(it_infected);

    // delete this from the infected vector
    infected.erase( it_infected );

    if (verbose)
        cout << "changing node status of infected " << this_infected << " to susceptible" << endl;

    // change node status of this node
    node_status[this_infected] = _S;

    if (verbose)
        cout << "erasing all SI edges which " << this_infected << " was part of" << endl;

    // erase all edges in the SI set which this infected is part of
    SI_edges.erase( 
            remove_if( 
                    SI_edges.begin(), 
                    SI_edges.end(),
                [&this_infected](const pair < size_t, size_t > & edge) { 
                        return edge.first == this_infected;
                }),
            SI_edges.end() 
        );

    if (verbose)
    {
        cout << "  SI edges after erasing = " << endl;
        print_SI_edges();
    }

    size_t const & this_susceptible = this_infected;

    if (verbose)
    {
        cout << "finding all edges which " << this_susceptible << " is part of and is connected to an infected " << endl;
        cout << "graph has size " << G->size() << endl;
        cout << "neighbor list to look at is "; 
        for(auto const &neighbor: (*G)[this_susceptible])
            cout << neighbor << " ";
        cout << endl;
    }

    // push the new SI edges
    for(auto const &neighbor: (*G)[this_susceptible])
    {
        if (verbose)
            cout << "new susceptible " << this_susceptible << " has neighbor " << neighbor << " with status " << node_status[neighbor] << endl;

        if (node_status[neighbor] == _I)
            SI_edges.push_back( make_pair(neighbor, this_susceptible) );
    }

    if (verbose)
    {
        cout << "  SI edges after pushing = " << endl;
        print_SI_edges();
    }

}

void Dyn_SIS::update_observables(
                double t
               )
{
    double _R0 = infection_rate * mean_degree / recovery_rate;
    R0.push_back(_R0);

    // compute SI
    SI.push_back(SI_edges.size());

    // compute I
    I.push_back(infected.size());

    // push back time
    time.push_back(t);
}

