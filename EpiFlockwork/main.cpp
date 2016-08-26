#include "SIS.h"

using namespace std;

int main()
{
    size_t N = 500;
    vector < tuple < size_t, size_t  > > E;
    double Q = 0.1;
    size_t t_run_total = 10000;
    double recovery_rate = 0.01;
    double rewiring_rate = 1.;
    size_t number_of_vaccinated = 0;
    size_t number_of_infected = 10;
    size_t seed = 1;
    double R0 = 1.;
    double k = 1. / (1.-Q);

    double infection_rate = R0 * recovery_rate / k;

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


    vector < pair < double, size_t > > I_of_t = get<0>(result);

    for(auto p: I_of_t)
    {
        cout << p.first << " " << p.second << endl;
    }

    //for(size_t i=0; i<result.size(); i++)
    //{
    //    cout << result[i]->at(0) << " " << result[i]->at(1) << "\n";
    //}


    return 0;
}
