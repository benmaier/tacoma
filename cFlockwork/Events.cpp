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
const size_t R = 2;
const size_t V = 3;

using namespace std;

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
		//cout << "node status: " << node_status[i] << " " << node_status[j] <<endl;
		//cout << "link: " << i << " " << j <<endl;
        throw domain_error( "There was a non SI link in the set of SI links. This should not happen." );
    }

    //change status of that node
    node_status[new_infected] = I;
    infected.insert(new_infected);

    //loop through the neighbors of i
    for(auto neigh_of_new_I : *G[new_infected] )
    {
        //get pair
        pair <size_t,size_t> current_edge = get_sorted_pair(new_infected,neigh_of_new_I);

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

    set< size_t >::const_iterator new_recovered_it(infected.begin());
    advance(new_recovered_it,element);

    size_t new_recovered = *new_recovered_it;
    //change status of that node
    node_status[new_recovered] = S;
    infected.erase(new_recovered);

    //loop through the neighbors of i
    for(auto neigh_of_new_S : *G[new_recovered] )
    {
        //get pair
        pair <size_t,size_t> current_edge = get_sorted_pair(new_recovered,neigh_of_new_S);

        //if the neighbor is infected, then the edge belongs in SI_E
        if (node_status[neigh_of_new_S] == I)
            SI_E.insert(current_edge);
        //if the neighbor is susceptible, then it does not belong in SI_E anymore
        else if (node_status[neigh_of_new_S] == S)
            SI_E.erase(current_edge);
    }

}

void SIR_recover(
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

    set< size_t >::const_iterator new_recovered_it(infected.begin());
    advance(new_recovered_it,element);

    size_t new_recovered = *new_recovered_it;
    //change status of that node
    node_status[new_recovered] = R;
    infected.erase(new_recovered);

    //loop through the neighbors of i
    for(auto neigh_of_new_R : *G[new_recovered] )
    {
        //if the neighbor is susceptible, then it doesnt belong in SI_E anymore
        if (node_status[neigh_of_new_R] == S)
        {
            //get edge
            pair <size_t,size_t> current_edge = get_sorted_pair(new_recovered,neigh_of_new_R);
            SI_E.erase(current_edge);
        }
    }
}

void SIRS_recover(
                 vector < set < size_t > * > & G, //Adjacency matrix
                 default_random_engine & generator, 
                 uniform_real_distribution<double> & distribution,
                 set < pair < size_t, size_t > > & SI_E, //edge list of SI links
                 vector < size_t > & node_status,
                 set < size_t > & infected,
                 set < size_t > & recovered
           )
{
    //get element of infected that recovers
    double r = distribution(generator);
    size_t element = infected.size() * r;

    set< size_t >::const_iterator new_recovered_it(infected.begin());
    advance(new_recovered_it,element);

    size_t new_recovered = *new_recovered_it;
    //change status of that node
    node_status[new_recovered] = R;
    infected.erase(new_recovered);
    recovered.insert(new_recovered);

    //loop through the neighbors of i
    for(auto neigh_of_new_R : *G[new_recovered] )
    {
        //if the neighbor is susceptible, then it doesnt belong in SI_E anymore
        //if (new_recovered==752) 
            //cout << "I->R, status 752 = " << node_status[new_recovered] << "; node_status[659] = "<< node_status[659] << endl;
        if (node_status[neigh_of_new_R] == S)
        {
            //get edge
            pair <size_t,size_t> current_edge = get_sorted_pair(new_recovered,neigh_of_new_R);
            SI_E.erase(current_edge);
        }
    }
}

void become_susceptible(
                 vector < set < size_t > * > & G, //Adjacency matrix
                 default_random_engine & generator, 
                 uniform_real_distribution<double> & distribution,
                 set < pair < size_t, size_t > > & SI_E, //edge list of SI links
                 vector < size_t > & node_status,
                 set < size_t > & recovered
           )
{
    //get element of infected that recovers
    double r = distribution(generator);
    size_t element = recovered.size() * r;

    set< size_t >::const_iterator new_susceptible_it(recovered.begin());
    advance(new_susceptible_it,element);

    size_t new_susceptible = *new_susceptible_it;
    //change status of that node
    node_status[new_susceptible] = S;
    recovered.erase(new_susceptible);

	//cout << "status neigh of new S" << endl;
    //loop through the neighbors of i
    for(auto neigh_of_new_S : *G[new_susceptible] )
    {
        //if the neighbor is infected, then the edge belongs in SI_E
		//cout << node_status[neigh_of_new_S] << endl;
        //if (new_susceptible==752) 
        //    cout << "status 752 = " << node_status[new_susceptible] << "; node_status[659] = "<< node_status[659] << endl;
        if (node_status[neigh_of_new_S] == I)
        {
            //get edge
            pair <size_t,size_t> current_edge = get_sorted_pair(new_susceptible,neigh_of_new_S);
            SI_E.insert(current_edge);
			//cout << "new_recovered: " << new_susceptible << " ; infected" << neigh_of_new_S <<endl;
			//cout << node_status[new_susceptible] << " " << node_status[neigh_of_new_S] << endl;
        }
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
    rewire_P(
                 vector < set < size_t > * > & G, //Adjacency matrix
                 double P,       //probability to connect with neighbors of neighbor
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
    size_t number_of_new_edges = 0;

    vector < pair < size_t, size_t > > edges_out;
    vector < pair < size_t, size_t > > edges_in;

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

    if ( distribution(generator) < P )
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
                 default_random_engine & generator, 
                 uniform_real_distribution<double> & distribution,
                 double & mean_degree,
                 set < pair < size_t, size_t > > & SI_E, //edge list of SI links
                 const vector < size_t > & node_status,
                 const bool use_preferential_node_selection
            )
{
    //choose two nodes
    size_t N = G.size();

    size_t i, j;

    if (use_preferential_node_selection)
    {
        i = arg_choose_from_vector(total_affinity, generator, distribution);
    }
    else
    {
        double r1 = distribution(generator);
        i = r1 * N;
    }

    size_t number_of_old_edges = 0;
    size_t number_of_new_edges = 0;

    vector < pair < size_t, size_t > > edges_out;
    vector < pair < size_t, size_t > > edges_in;

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

    if ( distribution(generator) < P )
    {
        // get node from social network structure
        size_t neighbor_index = arg_choose_from_vector(neighbor_affinity[i].second, generator, distribution);        
        j = neighbor_affinity[i].first[neighbor_index];

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

void random_rewire(
                 vector < set < size_t > * > & G, //Adjacency matrix
                 default_random_engine & generator, 
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

