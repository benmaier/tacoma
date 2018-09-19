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

#include "SIS_node_based.h"

using namespace std;
//namespace node_based_SIS = node_based_SIS;

void node_based_SIS::update_network(
    vector<set<size_t>> &_G,
    double t)
{
    // map this network to the current network
    G = &_G;

    // compute degree and R0
    size_t number_of_edges_times_two = 0;

    for (auto const &neighbors : *G)
        number_of_edges_times_two += neighbors.size();

    mean_degree = number_of_edges_times_two / (double)N;

    // update infected
    SI_count = 0;
    auto it_SI_degree = SI_degree.begin();
    for (auto &neighs : SI_Graph)
    {
        neighs.clear();
        *it_SI_degree = 0;
        ++it_SI_degree;
    }

    // for each infected, check its neighbors.
    // if the neighbor is susceptible, push it to the
    // endangered susceptible vector
    for (auto const &inf : infected)
        for (auto const &neighbor_of_infected : (*G)[inf])
            if (node_status[neighbor_of_infected] == EPI::S)
            {
                SI_Graph[inf].insert(neighbor_of_infected);
                ++SI_degree[inf];
                ++SI_count;
            }

    // update the arrays containing the observables
    update_observables(t);
}

void node_based_SIS::update_network(
    vector<set<size_t>> &_G,
    vector<pair<size_t, size_t>> &edges_in,
    vector<pair<size_t, size_t>> &edges_out,
    double t)
{
    // map this network to the current network
    G = &_G;

    // compute degree and R0
    size_t number_of_edges_times_two = 0;

    for (auto const &neighbors : *G)
        number_of_edges_times_two += neighbors.size();

    mean_degree = number_of_edges_times_two / (double)N;

    for (auto const &e : edges_out)
    {
        size_t u = e.first;
        size_t v = e.second;

        if ((node_status[u] == EPI::S) and (node_status[v] == EPI::I))
        {
            SI_Graph[v].erase(u);
            --SI_degree[v];
            --SI_count;
        }
        else if ((node_status[v] == EPI::S) and (node_status[u] == EPI::I))
        {
            SI_Graph[u].erase(v);
            --SI_degree[u];
            --SI_count;
        }
    }

    for (auto const &e : edges_in)
    {
        size_t u = e.first;
        size_t v = e.second;

        if ((node_status[u] == EPI::S) and (node_status[v] == EPI::I))
        {
            SI_Graph[v].insert(u);
            ++SI_degree[v];
            ++SI_count;
        }
        else if ((node_status[v] == EPI::S) and (node_status[u] == EPI::I))
        {
            SI_Graph[u].insert(v);
            ++SI_degree[u];
            ++SI_count;
        }
    }

    // update the arrays containing the observables
    update_observables(t);
}

void node_based_SIS::get_rates_and_Lambda(
    vector<double> &_rates,
    double &_Lambda)
{
    // delete current rates
    rates.clear();

    // compute rates of infection
    rates.push_back(infection_rate * SI_count);

    // compute rates of recovery
    rates.push_back(recovery_rate * infected.size());

    // return those new rates
    _rates = rates;
    _Lambda = accumulate(rates.begin(), rates.end(), 0.0);
}

void node_based_SIS::make_event(
    size_t const &event,
    double t)
{
    if (event == 0)
        infection_event();
    else if (event == 1)
        recovery_event();
    else
        throw length_error("node_based_SIS: chose event larger than rate vector which should not happen.");

    update_observables(t);
}

void node_based_SIS::infection_event()
{
    // choose random infected
    uniform_real_distribution<double> random_uni(0.0, 1.0);
    size_t this_infected = arg_choose_from_vector(SI_degree, generator, random_uni);

    // choose random susceptible neighbor of this infected
    uniform_int_distribution<size_t> random_susceptible(0, SI_Graph[this_infected].size() - 1);
    size_t this_susceptible_index = random_susceptible(generator);
    auto it_susceptible = SI_Graph[this_infected].begin();
    advance(it_susceptible, this_susceptible_index);
    size_t this_susceptible = *it_susceptible;

    // save this node as an infected
    infected.push_back(this_susceptible);

    // change node status of this node
    node_status[this_susceptible] = EPI::I;

    // erase all edges in the SI set where this susceptible is part of
    SI_Graph[this_infected].erase(this_susceptible);

    for (auto const &neighbor : (*G)[this_susceptible])
    {
        if (node_status[neighbor] == EPI::S)
        {
            SI_Graph[this_susceptible].insert(neighbor);
            ++SI_count;
            ++SI_degree[this_susceptible];
        }
        else if (node_status[neighbor] == EPI::I)
        {
            SI_Graph[neighbor].erase(this_susceptible);
            --SI_count;
            --SI_degree[neighbor];
        }
    }
}

void node_based_SIS::recovery_event()
{
    if ((prevent_disease_extinction) and (infected.size() == 1))
        return;

    // initialize uniform integer random distribution
    uniform_int_distribution<size_t> random_infected(0, infected.size() - 1);

    // find the index of the susceptible which will become infected
    size_t this_infected_index = random_infected(generator);
    auto it_infected = infected.begin() + this_infected_index;

    // get the node id of this infected about to be recovered
    size_t this_infected = *(it_infected);

    // delete this from the infected vector
    infected.erase(it_infected);

    // change node status of this node
    node_status[this_infected] = EPI::S;

    // erase all edges in the SI set which this infected is part of
    SI_Graph[this_infected].clear();
    SI_count -= SI_degree[this_infected];
    SI_degree[this_infected] = 0;

    for (auto &neigh : (*G)[this_infected])
        if (node_status[neigh] == EPI::I)
        {
            SI_Graph[neigh].insert(this_infected);
            ++SI_degree[neigh];
            ++SI_count;
        }
}

void node_based_SIS::update_observables(
    double t)
{
    double _R0 = infection_rate * mean_degree / recovery_rate;
    R0.push_back(_R0);

    // compute SI
    SI.push_back(SI_count);

    // compute I
    I.push_back(infected.size());

    // push back time
    time.push_back(t);
}
