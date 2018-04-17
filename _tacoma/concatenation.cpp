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
#include "conversion.h"
#include "social_trajectories.h"

using namespace std;

edge_changes
     concatenate_edge_changes(
            vector < edge_changes > & ec_lists,
            const bool verbose
        )
{
    edge_changes new_ec;

    auto ec_it = ec_lists.begin();

    // create graph
    size_t N = ec_it->N;

    bool all_equal_N = true;
    bool all_have_node_map = true;

    // reject the concatenation if N is different
    // in the different temporal networks
    // and there's no node mapping
    for(; ec_it != ec_lists.end(); ++ec_it)
    {
        all_equal_N = all_equal_N and (ec_it->N == N);
        all_have_node_map = all_have_node_map and (ec_it->int_to_node.size() == ec_it->N);
    }

    if ((not all_equal_N) and (not all_have_node_map))
        throw domain_error("Cannot reliably concatenate temporal networks with different node sizes where no sufficient node remaps are given.");

    bool remap_nodes = all_have_node_map;
    map < string, size_t > all_nodes_to_int;
    map < size_t, string > all_int_to_node;

    // construct the new node maps
    for(ec_it = ec_lists.begin(); ec_it != ec_lists.end(); ++ec_it)
    {
        for(auto const & this_i2n: ec_it->int_to_node)
        {
            string const &s = this_i2n.second;

            if (all_nodes_to_int.find(s) == all_nodes_to_int.end())
            {
                size_t new_int = all_nodes_to_int.size();
                all_nodes_to_int[s] = all_nodes_to_int.size();
                all_int_to_node[new_int] = s;
            }
        }
    }

    if (remap_nodes)
        N = all_nodes_to_int.size();


    // rewind
    ec_it = ec_lists.begin();
    double t0 = ec_it->t0;

    new_ec.N = N;
    new_ec.t0 = t0;
    new_ec.edges_initial = ec_it->edges_initial;


    vector < vector < pair < size_t, size_t > > > new_edges_in;
    vector < vector < pair < size_t, size_t > > > new_edges_out;

    vector < double > new_time;
    
    //vector < pair < size_t, size_t > > these_edges = ec_it->edges_initial;

    double last_tmax = t0;

    set < string > time_units;
    string all_notes = "";

    // TODO: write a procedure that efficiently checks the 
    // node integer maps and remaps integers which 
    // point to different node names
    // map < size_t, string > new_int_to_node;

    // iterate through all edge_changes
    while ( ec_it != ec_lists.end() )
    {

        time_units.insert( ec_it->time_unit );

        if (time_units.size() > 1)
            throw domain_error("instances of `edge_changes` have varying time_units");

        all_notes += "<< NEXT NOTES >> " + ec_it->notes + " ";


                
        vector < set < size_t > > G(N);
        vector < pair < size_t, size_t > > these_edges_initial((ec_it->edges_initial).begin(), (ec_it->edges_initial).begin()); 

        if (remap_nodes)
        {
            for(auto &edge: these_edges_initial)
            {
                size_t &i = edge.first;
                size_t &j = edge.second;

                i = all_nodes_to_int[ ec_it->int_to_node[i] ];
                j = all_nodes_to_int[ ec_it->int_to_node[j] ];

                if (i>j)
                    swap(i,j);
            }
        }
        graph_from_edgelist(G,these_edges_initial);

        // get time of current state

        // create iterators
        auto it_edges_in = (ec_it->edges_in).begin();
        auto it_edges_out = (ec_it->edges_out).begin();
        auto it_time = (ec_it->t).begin();

        while (     it_edges_in != (ec_it->edges_in).end() 
                and it_time != (ec_it->t).end() 
                and it_edges_out != (ec_it->edges_out).end() 
              )
        {

            // get a sorted edge list
            vector < pair < size_t, size_t > > these_edges_in(it_edges_in->begin(), it_edges_in->end());
            vector < pair < size_t, size_t > > these_edges_out(it_edges_out->begin(), it_edges_out->end());

            for(auto &edge: these_edges_out)
            {
                size_t &i = edge.first;
                size_t &j = edge.second;

                if (remap_nodes)
                {
                    i = all_nodes_to_int[ ec_it->int_to_node[i] ];
                    j = all_nodes_to_int[ ec_it->int_to_node[j] ];

                    if (i>j)
                        swap(i,j);
                }

                G[i].erase(j);
                G[j].erase(i);
            }

            for(auto &edge: these_edges_in)
            {
                size_t &i = edge.first;
                size_t &j = edge.second;

                if (remap_nodes)
                {
                    i = all_nodes_to_int[ ec_it->int_to_node[i] ];
                    j = all_nodes_to_int[ ec_it->int_to_node[j] ];

                    if (i>j)
                        swap(i,j);
                }

                G[i].insert(j);
                G[j].insert(i);
            }


            new_time.push_back( *it_time - (ec_it->t0) + last_tmax );
            new_edges_in.push_back(these_edges_in);
            new_edges_out.push_back(these_edges_out);
            
            // advance iterators
            it_edges_in++;
            it_edges_out++;
            it_time++;
        }

        // calculate edge changes between last state and next initial edge list
        if (ec_it + 1 != ec_lists.end())
        {
            vector < pair < size_t, size_t > > this_edge_list;
            edgelist_from_graph(this_edge_list,G);

            vector < pair < size_t, size_t > > next_edge_list(
                                                                ((ec_it+1)->edges_initial).begin(),
                                                                ((ec_it+1)->edges_initial).end()
                                                             );

            if (remap_nodes)
            {
                for(auto &edge: next_edge_list)
                {
                    size_t &i = edge.first;
                    size_t &j = edge.second;

                    i = all_nodes_to_int[ (ec_it+1)->int_to_node[i] ];
                    j = all_nodes_to_int[ (ec_it+1)->int_to_node[j] ];

                    if (i>j)
                        swap(i,j);
                }
            }

            // get a sorted edge list
            set < size_t > last_edge_integers = get_edge_integer_set(this_edge_list, N);
            set < size_t > edge_integers = get_edge_integer_set(next_edge_list, N);

            vector < size_t > incoming_edge_integers;
            vector < size_t > outgoing_edge_integers;

            // incoming edges are from the set (edges - last_edges) 
            set_difference(
                           edge_integers.begin(), edge_integers.end(),
                           last_edge_integers.begin(), last_edge_integers.end(),
                           back_inserter(incoming_edge_integers)
                          );

            // outgoing edges are from the set (last_edges - edges) 
            set_difference(
                           last_edge_integers.begin(), last_edge_integers.end(),
                           edge_integers.begin(), edge_integers.end(),
                           back_inserter(outgoing_edge_integers)
                          );

            vector < pair < size_t, size_t > > these_edges_in;
            vector < pair < size_t, size_t > > these_edges_out;

            for(auto const &e: incoming_edge_integers)
            {
                size_t j = e % N;
                size_t i = e / N;
                these_edges_in.push_back( make_pair(i, j) );
            }

            for(auto const &e: outgoing_edge_integers)
            {
                size_t j = e % N;
                size_t i = e / N;
                these_edges_out.push_back( make_pair(i, j) );
            }

            new_edges_in.push_back( these_edges_in );
            new_edges_out.push_back( these_edges_out );
            new_time.push_back(ec_it->tmax - ec_it->t0 + last_tmax );
        }

        last_tmax += ec_it->tmax - ec_it->t0;
        ec_it++;
    }

    double tmax = last_tmax;

    new_ec.tmax = tmax;
    new_ec.edges_in = new_edges_in;
    new_ec.edges_out = new_edges_out;
    new_ec.t = new_time;
    new_ec.notes = all_notes;
    new_ec.time_unit = *time_units.begin();
    new_ec.int_to_node = all_int_to_node;

    return new_ec;
}

edge_lists
     concatenate_edge_lists(
            vector < edge_lists > & el_lists,
            const bool verbose
            )
{
    edge_lists new_ec;

    auto ec_it = el_lists.begin();

    // create graph
    size_t N = ec_it->N;

    bool all_equal_N = true;
    bool all_have_node_map = true;

    // reject the concatenation if N is different
    // in the different temporal networks
    // and there's no node mapping
    for(; ec_it != el_lists.end(); ++ec_it)
    {
        all_equal_N = all_equal_N and (ec_it->N == N);
        all_have_node_map = all_have_node_map and (ec_it->int_to_node.size() == ec_it->N);
    }

    if ((not all_equal_N) and (not all_have_node_map))
        throw domain_error("Cannot reliably concatenate temporal networks with different node sizes where no sufficient node remaps are given.");

    bool remap_nodes = all_have_node_map;
    map < string, size_t > all_nodes_to_int;
    map < size_t, string > all_int_to_node;

    // construct the new node maps
    for(ec_it = el_lists.begin(); ec_it != el_lists.end(); ++ec_it)
    {
        for(auto const & this_i2n: ec_it->int_to_node)
        {
            string const &s = this_i2n.second;

            if (all_nodes_to_int.find(s) == all_nodes_to_int.end())
            {
                size_t new_int = all_nodes_to_int.size();
                all_nodes_to_int[s] = all_nodes_to_int.size();
                all_int_to_node[new_int] = s;
            }
        }
    }

    if (remap_nodes)
        N = all_nodes_to_int.size();

    // rewind
    ec_it = el_lists.begin();

    vector < vector < pair < size_t, size_t > > > new_edges;
    vector < double > new_time;
    
    double last_tmax = ec_it->t.front();

    set < string > time_units;
    string all_notes = "";

    // TODO: write a procedure that efficiently checks the 
    // node integer maps and remaps integers which 
    // point to different node names

    // iterate through all edge_changes
    while ( ec_it != el_lists.end() )
    {
        time_units.insert( ec_it->time_unit );

        if (time_units.size() > 1)
            throw domain_error("instances of `edge_lists` have varying time_units");

        all_notes += "<< NEXT NOTES >> " + ec_it->notes + " ";
                
        // create iterators
        auto it_edges = (ec_it->edges).begin();
        auto it_time = (ec_it->t).begin();

        while (     
                    it_edges != (ec_it->edges).end() 
                and it_time != (ec_it->t).end() 
              )
        {
            vector < pair < size_t, size_t > > these_edges( it_edges->begin(), it_edges->end() );

            if (remap_nodes)
            {
                for(auto &edge: these_edges)
                {
                    size_t &i = edge.first;
                    size_t &j = edge.second;

                    i = all_nodes_to_int[ ec_it->int_to_node[i] ];
                    j = all_nodes_to_int[ ec_it->int_to_node[j] ];

                    if (i>j)
                        swap(i,j);
                }
            }

            new_edges.push_back(these_edges);
            new_time.push_back( *it_time - (ec_it->t.front()) + last_tmax );

            // advance iterators
            it_edges++;
            it_time++;
        }

        last_tmax += ec_it->tmax - ec_it->t.front();
        ec_it++;
    }

    double tmax = last_tmax;

    new_ec.N = N;
    new_ec.tmax = tmax;
    new_ec.edges = new_edges;
    new_ec.t = new_time;
    new_ec.notes = all_notes;
    new_ec.time_unit = *time_units.begin();
    new_ec.int_to_node = all_int_to_node;


    return new_ec;
}
