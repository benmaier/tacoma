#include "Events.h"
#include "Utilities.h"
#include <ctime>

using namespace std;

int main()
{
    vector < set < size_t > * > G;
    vector < size_t > node_ints;
    vector < size_t > node_status;
    set < tuple<size_t,size_t> > SI_E;

    size_t N = 10;

    for(size_t node=0; node<N; node++)
    {
        G.push_back(new set <size_t>);
        node_ints.push_back(node);
        node_status.push_back(0);
    }

    for(size_t i=0; i<N-1; i++)
        for(size_t j=i+1; j<N; j++)
        {
            G[i]->insert(j);
            G[j]->insert(i);
        }


    //initialize random generators
    default_random_engine generator(time(nullptr));
    uniform_real_distribution<double> uni_distribution(0.,1.);

    random_rewire(G,generator,uni_distribution,SI_E,node_status,node_ints);
    random_rewire(G,generator,uni_distribution,SI_E,node_status,node_ints);
    random_rewire(G,generator,uni_distribution,SI_E,node_status,node_ints);
    random_rewire(G,generator,uni_distribution,SI_E,node_status,node_ints);
    random_rewire(G,generator,uni_distribution,SI_E,node_status,node_ints);
    //choose_random_unique(node_ints.begin(),node_ints.end(),3,generator,uni_distribution);

    for(size_t i=0; i<N; i++)
        for(auto neigh: *G[i])
        {
            cout << i << " " << neigh << endl;
        }

    for(size_t i=0; i<N; i++)
    {
        cout << node_ints[i] << endl;
    }

    return 0;
}
