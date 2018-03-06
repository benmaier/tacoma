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

#include "Utilities.h"
#include "ResultClasses.h"
#include "dyn_RGG.h"

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

using namespace std;

void dyn_RGG_update_positions(vector < pair < double, double > > &pos,
                      default_random_engine &generator,
                      uniform_real_distribution<double> &uni_distribution,
                      const double &step_distance)
{
    const size_t N = pos.size();
    for(size_t node = 0; node < N; node++)
    {
        double theta = uni_distribution(generator)*2*M_PI;
        pos[node].first += step_distance * cos(theta);
        pos[node].second += step_distance * sin(theta);

        if (pos[node].first < 0)
            pos[node].first += 1.0;
        if (pos[node].first >= 1.0)
            pos[node].first -= 1.0;
        if (pos[node].second < 0)
            pos[node].second += 1.0;
        if (pos[node].second >= 1.0)
            pos[node].second -= 1.0;
    }
}

void dyn_RGG_write_edge_list(
                   vector < pair < size_t, size_t > > & edges,
                   vector < pair < double, double > > const &pos,
                   const double &R,
                   const bool &PBC_distance
        )
{
    size_t const N = pos.size();
    int copymin, copymax;
    if (PBC_distance)
    {
        copymin = -1;
        copymax = +1;
    }
    else
    {
        copymin = 0;
        copymax = 0;
    }
    for(size_t i=0; i<N-1; ++i)
    {
        for(size_t j=i+1; j<N; ++j)
        {
            for(int ix = copymin; ix <= copymax; ix++)
            {
                for(int iy = copymin; iy <= copymax; iy++)
                {                    
                    double const x1 = pos[i].first;
                    double const y1 = pos[i].second;
                    double const x2 = pos[j].first + ix;
                    double const y2 = pos[j].second + iy;
                    double distance = sqrt( pow(x1-x2,2) + pow(y1-y2,2));
                    if (distance < R)
                    {
                        edges.push_back(make_pair(i,j));
                    }
                }
            }
        }
    }
}

edge_lists_with_histograms
     dynamic_RGG(
             const size_t N,       //number of nodes
             const size_t t_run_total,
             double step_distance,
             const double mean_link_duration,
             const double critical_density,
             const bool   PBC_distance,
             const bool   record_sizes_and_durations,
             const size_t seed,
             const bool verbose
        )
{

    const double R = sqrt( critical_density / (double) N );

    if ( (step_distance<=0.0) && (mean_link_duration<=0.0) )
        throw domain_error("Please provide either positive step_distance xor positive mean_link_duration.");
    else if ( (step_distance>0.0) && (mean_link_duration>0.0) )
        throw domain_error("Please provide either positive step_distance xor positive mean_link_duration, not both positive.");
    else if (mean_link_duration>0.0)
        step_distance = 1.12 * R / ( mean_link_duration - 0.26 );

    //initialize random generators
    default_random_engine generator;
    if (seed == 0)
        randomly_seed_engine(generator);
    else
        generator.seed(seed);

    uniform_real_distribution<double> uni_distribution(0.,1.);

    vector < pair < double, double > > pos;

    //intialize positions
    for(size_t node=0; node<N; node++)
    {
        double x = uni_distribution(generator);
        double y = uni_distribution(generator);
        pos.push_back(make_pair(x,y));
    }

    if (verbose)
    {
        cout << "initialized" << endl;
        cout << "R = " << R << endl;
        cout << "step_distance = " << step_distance << endl;
    }


    vector < vector < pair <size_t,size_t> > > edges;
    vector < double > time;

    // maps and sets for measuring edge durations
    set < size_t > initial_edges;
    map < size_t, size_t > current_edges;
    map < size_t, size_t >::iterator current_edge_iterator;
    vector < size_t > group_durations;
    vector < map < size_t, size_t > > histograms;



    for(size_t t = 0; t<t_run_total; t++)
    {
        if (verbose)
        {
            cout << " ============== " << endl;
            cout << "t = " << t << endl;
        }

        vector < pair <size_t,size_t> > these_edges;

        dyn_RGG_update_positions(pos, generator, uni_distribution, step_distance);

        dyn_RGG_write_edge_list(these_edges, pos, R, PBC_distance);

        edges.push_back(these_edges);

        if (record_sizes_and_durations)
        {
            // check for new edges or terminated edges
            if (t == 0)
            {
                for(auto const &edge: these_edges)
                {
                    size_t edge_int = get_edge_int(edge,N);
                    initial_edges.insert(edge_int);
                    current_edges[edge_int] = t;
                }
            }
            else
            {
                // check if any of the new edges are actually new
                // by checking wether or not they are in the current_edges dictionary
                for(auto const &edge: these_edges)
                {
                    size_t edge_int = get_edge_int(edge,N);
                    const bool already_in = current_edges.find(edge_int) != current_edges.end();
                    if (not already_in)
                    {
                        current_edges[edge_int] = t;
                    }
                }

                // check if any of the old edges are not anymore in the new edges
                vector < size_t > edges_to_delete; 
                for(current_edge_iterator = current_edges.begin();
                    current_edge_iterator != current_edges.end(); 
                    current_edge_iterator++)
                {
                    size_t edge_int = current_edge_iterator->first;
                    size_t i = edge_int / N;
                    size_t j = edge_int % N;
                    pair < size_t, size_t > edge = make_pair(i,j);

                    const bool not_in = find(these_edges.begin(),
                                             these_edges.end(),
                                             edge
                                             ) == these_edges.end();
                    if (not_in)
                    {
                        const bool is_initial_edge =    initial_edges.find(edge_int) 
                                                     != initial_edges.end();
                        if (initial_edges.size()==0 ||
                            not is_initial_edge
                           )
                        {
                            size_t duration = t - current_edge_iterator->second;
                            group_durations.push_back(duration);
                        }
                        else if (is_initial_edge) {
                            initial_edges.erase(edge_int);
                        }
                        edges_to_delete.push_back(edge_int);
                    }
                }
                for( auto const &edge_int: edges_to_delete)
                {
                    current_edges.erase(edge_int);
                }
            }

            // get current
            map < size_t, size_t > counter;

            get_component_size_histogram_from_edgelist(N, these_edges, counter);
            histograms.push_back(counter);
        }

        time.push_back(t);
    }

    edge_lists_with_histograms result;

    result.t = time;
    result.N = N;
    result.tmax = t_run_total;
    result.edges = edges;
    result.size_histograms = histograms;
    result.group_durations = group_durations;

    return result;
}

