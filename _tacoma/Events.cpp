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
#include "Events.h"

const size_t S = 0;
const size_t I = 1;
const bool EVENT_VERBOSE = false;

using namespace std;

void rewire(
                 vector < set < size_t > * > & G, //Adjacency matrix
                 double Q,       //probability to connect with neighbors of neighbor
                 mt19937_64 & generator, 
                 uniform_real_distribution<double> & distribution,
                 double & mean_degree,
                 set < pair < size_t, size_t > > & SI_E, //edge list of SI links
                 const vector < size_t > & node_status
            )
{
    //choose two nodes
    size_t N = G.size();
    size_t i,j;

    uniform_int_distribution<size_t> node1(0,N-1);
    uniform_int_distribution<size_t> node2(0,N-2);

    i = node1(generator);
    j = node2(generator);

    if (j>=i)
        j++;


    size_t number_of_old_edges = 0;
    size_t number_of_new_edges = 0;

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
            SI_E.erase( get_sorted_pair(i,neigh_i) );
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
                SI_E.insert( get_sorted_pair(i,neigh_j) );
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
        SI_E.insert( get_sorted_pair(i,j) );
    }

    //calculate new number of edges
    mean_degree += 2.0 * ( double(number_of_new_edges) - double(number_of_old_edges) ) / double(N);

}

pair < vector < pair < size_t, size_t > >, vector < pair < size_t, size_t > > > 
    rewire_P_without_SI_checking(
                 vector < set < size_t > * > & G, //Adjacency matrix
                 double P,       //probability to connect with neighbors of neighbor
                 mt19937_64 & generator, 
                 uniform_real_distribution<double> & distribution
            )
{
    //choose two nodes
    size_t N = G.size();
    size_t i,j;

    uniform_int_distribution<size_t> node1(0,N-1);
    uniform_int_distribution<size_t> node2(0,N-2);

    i = node1(generator);
    j = node2(generator);

    if (j>=i)
        j++;

    bool do_rewiring = distribution(generator) < P;

    vector < pair < size_t, size_t > > edges_out;
    vector < pair < size_t, size_t > > edges_in;

    //check if new neighbor is actually an old neighbor
    //and if this is the case return an empty event
    if ( do_rewiring and (G[i]->find(j) != G[i]->end()) )
    {
        return make_pair(edges_out,edges_in);
    }

    //loop through the neighbors of i
    for(auto neigh_i : *G[i] )
    {
        //and erase the link to i
        G[neigh_i]->erase(i);
        
        pair < size_t, size_t > current_edge = get_sorted_pair(i,neigh_i);
        edges_out.push_back( current_edge );
    } 

    //erase the links from the perspective of i
    G[i]->clear();

    if ( do_rewiring )
    {
        //loop through the neighbors of j
        for(auto neigh_j : *G[j] ) 
        {
            G[ neigh_j ]->insert( i );
            G[ i ]->insert( neigh_j );

            pair < size_t, size_t > current_edge = get_sorted_pair(i,neigh_j);
            edges_in.push_back( current_edge );
        }

        //add j as neighbor and count additional edge
        G[ j ]->insert( i );
        G[ i ]->insert( j );

        pair < size_t, size_t > current_edge = get_sorted_pair(i,j);
        edges_in.push_back( current_edge );
    }

    return make_pair(edges_out,edges_in);
}

pair < vector < pair < size_t, size_t > >, vector < pair < size_t, size_t > > > 
    rewire_P_without_SI_checking(
                 vector < set < size_t > > & G, //Adjacency matrix
                 double P,       //probability to connect with neighbors of neighbor
                 mt19937_64 & generator, 
                 uniform_real_distribution<double> & distribution
            )
{
    //choose two nodes
    size_t N = G.size();
    size_t i,j;

    uniform_int_distribution<size_t> node1(0,N-1);
    uniform_int_distribution<size_t> node2(0,N-2);

    i = node1(generator);
    j = node2(generator);

    if (j>=i)
        j++;

    bool do_rewiring = distribution(generator) < P;

    vector < pair < size_t, size_t > > edges_out;
    vector < pair < size_t, size_t > > edges_in;

    //check if new neighbor is actually an old neighbor
    //and if this is the case return an empty event
    if ( do_rewiring and (G[i].find(j) != G[i].end()) )
    {
        return make_pair(edges_out,edges_in);
    }

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

    return make_pair(edges_out,edges_in);
}



pair < vector < pair < size_t, size_t > >, vector < pair < size_t, size_t > > > 
    rewire_P_without_SI_checking_single_node(
                 size_t i,
                 vector < set < size_t > * > & G, //Adjacency matrix
                 double P,       //probability to connect with neighbors of neighbor
                 mt19937_64 & generator, 
                 uniform_real_distribution<double> & distribution
            )
{
    //choose second node
    size_t N = G.size();

    uniform_int_distribution<size_t> random_node(0,N-2);
    size_t j = random_node(generator);

    if (j>=i)
        ++j;

    //cout << "found pair (" << i << ", " << j << ")" << endl;

    bool do_rewiring = distribution(generator) < P;

    vector < pair < size_t, size_t > > edges_out;
    vector < pair < size_t, size_t > > edges_in;

    //check if new neighbor is actually an old neighbor
    //and if this is the case return an empty event
    if ( do_rewiring and (G[i]->find(j) != G[i]->end()) )
    {
        return make_pair(edges_out,edges_in);
    }

    //loop through the neighbors of i
    for(auto neigh_i : *G[i] )
    {
        //and erase the link to i
        G[neigh_i]->erase(i);
        
        pair < size_t, size_t > current_edge = get_sorted_pair(i,neigh_i);
        edges_out.push_back( current_edge );
    } 

    //erase the links from the perspective of i
    G[i]->clear();

    if ( do_rewiring )
    {
        //loop through the neighbors of j
        for(auto neigh_j : *G[j] ) 
        {
            G[ neigh_j ]->insert( i );
            G[ i ]->insert( neigh_j );

            pair < size_t, size_t > current_edge = get_sorted_pair(i,neigh_j);
            edges_in.push_back( current_edge );
        }

        //add j as neighbor and count additional edge
        G[ j ]->insert( i );
        G[ i ]->insert( j );

        pair < size_t, size_t > current_edge = get_sorted_pair(i,j);
        edges_in.push_back( current_edge );
    }

    return make_pair(edges_out,edges_in);
}

pair < vector < pair < size_t, size_t > >, vector < pair < size_t, size_t > > > 
    rewire_P_without_SI_checking_single_node_neighbor_affinity(
                 size_t i,
                 vector < set < size_t > * > & G, //Adjacency matrix
                 double P,       //probability to connect with neighbors of neighbor
                 vector < pair < vector < size_t >, vector < double > > > &neighbor_affinity,
                 mt19937_64 & generator, 
                 uniform_real_distribution<double> & distribution
            )
{
    //choose second node
    size_t N = G.size();


    bool do_rewiring = neighbor_affinity[i].second.size() > 0 && (distribution(generator) < P);
    size_t j;

    if (do_rewiring)
    {
        size_t neighbor_index = arg_choose_from_vector(neighbor_affinity[i].second, generator, distribution);
        j = neighbor_affinity[i].first[neighbor_index];
    }

    //cout << "found pair (" << i << ", " << j << ")" << endl;


    vector < pair < size_t, size_t > > edges_out;
    vector < pair < size_t, size_t > > edges_in;

    //check if new neighbor is actually an old neighbor
    //and if this is the case return an empty event
    if ( do_rewiring and (G[i]->find(j) != G[i]->end()) )
    {
        return make_pair(edges_out,edges_in);
    }

    //loop through the neighbors of i
    for(auto neigh_i : *G[i] )
    {
        //and erase the link to i
        G[neigh_i]->erase(i);
        
        pair < size_t, size_t > current_edge = get_sorted_pair(i,neigh_i);
        edges_out.push_back( current_edge );
    } 

    //erase the links from the perspective of i
    G[i]->clear();

    if ( do_rewiring )
    {
        //loop through the neighbors of j
        for(auto neigh_j : *G[j] ) 
        {
            G[ neigh_j ]->insert( i );
            G[ i ]->insert( neigh_j );

            pair < size_t, size_t > current_edge = get_sorted_pair(i,neigh_j);
            edges_in.push_back( current_edge );
        }

        //add j as neighbor and count additional edge
        G[ j ]->insert( i );
        G[ i ]->insert( j );

        pair < size_t, size_t > current_edge = get_sorted_pair(i,j);
        edges_in.push_back( current_edge );
    }

    return make_pair(edges_out,edges_in);
}



pair < vector < pair < size_t, size_t > >, vector < pair < size_t, size_t > > > 
    rewire_P(
                 vector < set < size_t > * > & G, //Adjacency matrix
                 double P,       //probability to connect with neighbors of neighbor
                 mt19937_64 & generator, 
                 uniform_real_distribution<double> & distribution,
                 double & mean_degree,
                 set < pair < size_t, size_t > > & SI_E, //edge list of SI links
                 const vector < size_t > & node_status
            )
{
    //choose two nodes
    size_t N = G.size();
    size_t i,j;

    uniform_int_distribution<size_t> node1(0,N-1);
    uniform_int_distribution<size_t> node2(0,N-2);

    i = node1(generator);
    j = node2(generator);

    if (j>=i)
        j++;

    bool do_rewiring = distribution(generator) < P;

    size_t number_of_old_edges = 0;
    size_t number_of_new_edges = 0;

    vector < pair < size_t, size_t > > edges_out;
    vector < pair < size_t, size_t > > edges_in;

    //check if new neighbor is actually an old neighbor
    //and if this is the case return an empty event
    if ( do_rewiring and (G[i]->find(j) != G[i]->end()) )
    {
        return make_pair(edges_out,edges_in);
    }

    //loop through the neighbors of i
    for(auto neigh_i : *G[i] )
    {
        //and erase the link to i
        G[neigh_i]->erase(i);
        number_of_old_edges++;
        
        pair < size_t, size_t > current_edge = get_sorted_pair(i,neigh_i);
        edges_out.push_back( current_edge );

        //check if we erased an SI link
        if ( 
             ( (node_status[i] == I) && (node_status[neigh_i] == S) ) ||
             ( (node_status[i] == S) && (node_status[neigh_i] == I) )
           ) 
        {
            SI_E.erase( current_edge );
        }
    } 

    //erase the links from the perspective of i
    G[i]->clear();

    if ( do_rewiring )
    {
        //loop through the neighbors of j
        for(auto neigh_j : *G[j] ) 
        {
            G[ neigh_j ]->insert( i );
            G[ i ]->insert( neigh_j );
            number_of_new_edges++;

            pair < size_t, size_t > current_edge = get_sorted_pair(i,neigh_j);
            edges_in.push_back( current_edge );

            //check if we created an SI link
            if ( 
                 ( (node_status[i] == I) && (node_status[neigh_j] == S) ) ||
                 ( (node_status[i] == S) && (node_status[neigh_j] == I) )
               ) 
            {
                SI_E.insert( current_edge );
            }
        }

        //add j as neighbor and count additional edge
        G[ j ]->insert( i );
        G[ i ]->insert( j );
        number_of_new_edges++;

        pair < size_t, size_t > current_edge = get_sorted_pair(i,j);
        edges_in.push_back( current_edge );

        //check if we created an SI link
        if ( 
             ( (node_status[i] == I) && (node_status[j] == S) ) ||
             ( (node_status[i] == S) && (node_status[j] == I) )
           ) 
        {
            SI_E.insert( current_edge );
        }
    }

    //calculate new number of edges
    mean_degree += 2.0 * ( double(number_of_new_edges) - double(number_of_old_edges) ) / double(N);

    return make_pair(edges_out,edges_in);

}

pair < vector < pair < size_t, size_t > >, vector < pair < size_t, size_t > > > 
    rewire_P_neighbor_affinity(
                 vector < set < size_t > * > & G, //Adjacency matrix
                 double P,       //probability to connect with neighbors of neighbor
                 vector < pair < vector < size_t >, vector < double > > > &neighbor_affinity,
                 vector < double > &total_affinity,
                 mt19937_64 & generator, 
                 uniform_real_distribution<double> & distribution,
                 double & mean_degree,
                 set < pair < size_t, size_t > > & SI_E, //edge list of SI links
                 const vector < size_t > & node_status,
                 const bool use_preferential_node_selection
            )
{
    if (EVENT_VERBOSE)
        cout << "entering event function" << endl;
    //choose two nodes
    size_t N = G.size();

    size_t i, j;

    bool do_rewiring = distribution(generator) < P;


    if (EVENT_VERBOSE)
        cout << "choosing node" << endl;

    size_t number_of_old_edges = 0;
    size_t number_of_new_edges = 0;

    vector < pair < size_t, size_t > > edges_out;
    vector < pair < size_t, size_t > > edges_in;

    if (use_preferential_node_selection)
    {
        i = arg_choose_from_vector(total_affinity, generator, distribution);
    }
    else
    {   
        do
        {
            double r1 = distribution(generator);
            i = r1 * N;
        } while ( neighbor_affinity[i].second.size() == 0 );
        
    }

    // get second node from social network structure
    if (EVENT_VERBOSE)
        cout << "getting neighbor node of neighbor " << i << endl;


    size_t neighbor_index = arg_choose_from_vector(neighbor_affinity[i].second, generator, distribution);        

    if (EVENT_VERBOSE)
        cout << "found neighbor index... " << neighbor_index << endl;

    j = neighbor_affinity[i].first[neighbor_index];

    if (EVENT_VERBOSE)
        cout << "found neighbor node... " << j << endl;

    //check if new neighbor is actually an old neighbor
    //and if this is the case return an empty event
    if ( do_rewiring and (G[i]->find(j) != G[i]->end()) )
    {
        return make_pair(edges_out,edges_in);
    }

    if (EVENT_VERBOSE)
        cout << "chose node " << i << endl;

    if (EVENT_VERBOSE)
        cout << "looping through neighbors of " << i << endl;

    //loop through the neighbors of i
    for(auto neigh_i : *G[i] )
    {
        //and erase the link to i
        if (EVENT_VERBOSE)
            cout << "erasing " << i << " from list of neighbors of node " << neigh_i <<  endl;

        G[neigh_i]->erase(i);
        number_of_old_edges++;
        
        pair < size_t, size_t > current_edge = get_sorted_pair(i,neigh_i);
        edges_out.push_back( current_edge );

        //check if we erased an SI link
        if ( 
             ( (node_status[i] == I) && (node_status[neigh_i] == S) ) ||
             ( (node_status[i] == S) && (node_status[neigh_i] == I) )
           ) 
        {
            SI_E.erase( current_edge );
        }
    } 

    if (EVENT_VERBOSE)
        cout << "clearing neighbors of " << i << endl;

    //erase the links from the perspective of i
    G[i]->clear();

    if ( do_rewiring )
    {

        //loop through the neighbors of j
        for(auto neigh_j : *G[j] ) 
        {
            G[ neigh_j ]->insert( i );
            G[ i ]->insert( neigh_j );
            number_of_new_edges++;

            pair < size_t, size_t > current_edge = get_sorted_pair(i,neigh_j);
            edges_in.push_back( current_edge );

            //check if we created an SI link
            if ( 
                 ( (node_status[i] == I) && (node_status[neigh_j] == S) ) ||
                 ( (node_status[i] == S) && (node_status[neigh_j] == I) )
               ) 
            {
                SI_E.insert( current_edge );
            }
        }

        //add j as neighbor and count additional edge
        G[ j ]->insert( i );
        G[ i ]->insert( j );
        number_of_new_edges++;

        pair < size_t, size_t > current_edge = get_sorted_pair(i,j);
        edges_in.push_back( current_edge );

        //check if we created an SI link
        if ( 
             ( (node_status[i] == I) && (node_status[j] == S) ) ||
             ( (node_status[i] == S) && (node_status[j] == I) )
           ) 
        {
            SI_E.insert( current_edge );
        }
    }

    //calculate new number of edges
    mean_degree += 2.0 * ( double(number_of_new_edges) - double(number_of_old_edges) ) / double(N);

    if (EVENT_VERBOSE)
        cout << "leaving event function" << endl;

    return make_pair(edges_out,edges_in);

}

pair < vector < pair < size_t, size_t > >, vector < pair < size_t, size_t > > > 
    rewire_P_neighbor_affinity_without_SI_checking(
                 vector < set < size_t > * > & G, //Adjacency matrix
                 double P,       //probability to connect with neighbors of neighbor
                 vector < pair < vector < size_t >, vector < double > > > &neighbor_affinity,
                 vector < double > &total_affinity,
                 mt19937_64 & generator, 
                 uniform_real_distribution<double> & distribution,
                 const bool use_preferential_node_selection
            )
{
    if (EVENT_VERBOSE)
        cout << "entering event function" << endl;
    //choose two nodes
    size_t N = G.size();

    size_t i, j;

    bool do_rewiring = distribution(generator) < P;


    if (EVENT_VERBOSE)
        cout << "choosing node" << endl;

    vector < pair < size_t, size_t > > edges_out;
    vector < pair < size_t, size_t > > edges_in;

    if (use_preferential_node_selection)
    {
        i = arg_choose_from_vector(total_affinity, generator, distribution);
    }
    else
    {   
        do
        {
            double r1 = distribution(generator);
            i = r1 * N;
        } while ( neighbor_affinity[i].second.size() == 0 );
        
    }

    // get second node from social network structure
    if (EVENT_VERBOSE)
        cout << "getting neighbor node of neighbor " << i << endl;


    size_t neighbor_index = arg_choose_from_vector(neighbor_affinity[i].second, generator, distribution);        

    if (EVENT_VERBOSE)
        cout << "found neighbor index... " << neighbor_index << endl;

    j = neighbor_affinity[i].first[neighbor_index];

    if (EVENT_VERBOSE)
        cout << "found neighbor node... " << j << endl;

    //check if new neighbor is actually an old neighbor
    //and if this is the case return an empty event
    if ( do_rewiring and (G[i]->find(j) != G[i]->end()) )
    {
        return make_pair(edges_out,edges_in);
    }

    if (EVENT_VERBOSE)
        cout << "chose node " << i << endl;

    if (EVENT_VERBOSE)
        cout << "looping through neighbors of " << i << endl;

    //loop through the neighbors of i
    for(auto neigh_i : *G[i] )
    {
        //and erase the link to i
        if (EVENT_VERBOSE)
            cout << "erasing " << i << " from list of neighbors of node " << neigh_i <<  endl;

        G[neigh_i]->erase(i);
        
        pair < size_t, size_t > current_edge = get_sorted_pair(i,neigh_i);
        edges_out.push_back( current_edge );
    } 

    if (EVENT_VERBOSE)
        cout << "clearing neighbors of " << i << endl;

    //erase the links from the perspective of i
    G[i]->clear();

    if ( do_rewiring )
    {

        //loop through the neighbors of j
        for(auto neigh_j : *G[j] ) 
        {
            G[ neigh_j ]->insert( i );
            G[ i ]->insert( neigh_j );

            pair < size_t, size_t > current_edge = get_sorted_pair(i,neigh_j);
            edges_in.push_back( current_edge );

        }

        //add j as neighbor and count additional edge
        G[ j ]->insert( i );
        G[ i ]->insert( j );

        pair < size_t, size_t > current_edge = get_sorted_pair(i,j);
        edges_in.push_back( current_edge );

    }

    //calculate new number of edges
    if (EVENT_VERBOSE)
        cout << "leaving event function" << endl;

    return make_pair(edges_out,edges_in);

}

void random_rewire(
                 vector < set < size_t > * > & G, //Adjacency matrix
                 mt19937_64 & generator, 
                 uniform_real_distribution<double> & distribution,
                 set < pair < size_t, size_t > > & SI_E, //edge list of SI links
                 const vector < size_t > & node_status,
                 vector < size_t > & node_ints
            )
{
    //choose two nodes
    size_t N = G.size();
    size_t i = N * distribution(generator);

    size_t number_of_old_edges = 0;

    //loop through the neighbors of i
    for(auto neigh : *G[i] )
    {
        //and erase the link to i
        G[neigh]->erase(i);
        number_of_old_edges++;

        //check if we erased an SI link
        if ( 
             ( (node_status[i] == I) && (node_status[neigh] == S) ) ||
             ( (node_status[i] == S) && (node_status[neigh] == I) )
           ) 
        {
            SI_E.erase( get_sorted_pair(i,neigh) );
        }
    } 

    if (number_of_old_edges > 0)
    {
        //erase the links from the perspective of i
        G[i]->clear();

        choose_random_unique(node_ints.begin(), node_ints.end(), number_of_old_edges+1,generator,distribution);

        //loop through the new_neighbors
        size_t neigh_int = 0;
        while (neigh_int < number_of_old_edges)
        {
            size_t neigh = node_ints[neigh_int];

            //add as neighbor of i if it's not i itself
            if ( neigh != i )
            {
                G[ neigh ]->insert( i );
                G[ i ]->insert( neigh );

                //check if we created an SI link
                if ( 
                     ( (node_status[i] == I) && (node_status[neigh] == S) ) ||
                     ( (node_status[i] == S) && (node_status[neigh] == I) )
                   ) 
                {
                    SI_E.insert( get_sorted_pair(i,neigh) );
                }
            } else {
                //if it's i itself, skip this neighbor and increase the total
                //number of evaluated neighbors by one
                number_of_old_edges++; 
            }

            //advance to next neighbor
            neigh_int++;
        }
    }
}

void random_rewire_without_SI_checking(
                 vector < set < size_t > * > & G, //Adjacency matrix
                 mt19937_64 & generator, 
                 uniform_real_distribution<double> & distribution,
                 vector < size_t > & node_ints
            )
{
    //choose two nodes
    size_t N = G.size();
    size_t i = N * distribution(generator);

    size_t number_of_old_edges = 0;

    //loop through the neighbors of i
    for(auto neigh : *G[i] )
    {
        //and erase the link to i
        G[neigh]->erase(i);
        number_of_old_edges++;
    } 

    if (number_of_old_edges > 0)
    {
        //erase the links from the perspective of i
        G[i]->clear();

        choose_random_unique(node_ints.begin(), node_ints.end(), number_of_old_edges+1,generator,distribution);

        //loop through the new_neighbors
        size_t neigh_int = 0;
        while (neigh_int < number_of_old_edges)
        {
            size_t neigh = node_ints[neigh_int];

            //add as neighbor of i if it's not i itself
            if ( neigh != i )
            {
                G[ neigh ]->insert( i );
                G[ i ]->insert( neigh );

            } else {
                //if it's i itself, skip this neighbor and increase the total
                //number of evaluated neighbors by one
                number_of_old_edges++; 
            }

            //advance to next neighbor
            neigh_int++;
        }
    }
}

