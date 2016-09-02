#include "SIS.h"
#include "SIR.h"
#include "SIRS.h"
#include "ResultClasses.h"

using namespace std;

int main()
{
    size_t N = 1000;
    vector < pair < size_t, size_t  > > E;
    double Q = 0.9;
    size_t t_run_total = 1000;
    double recovery_rate = 1.0;
    double rewiring_rate = 1.0;
    double susceptible_rate = 1.0;
    size_t number_of_vaccinated = 10;
    size_t number_of_infected = 10;
    size_t seed = 1;
    double R0 = 1.5;
    double k = 1. / (1.-Q);
    double p = k / (N-1);
    bool use_random_rewiring = false;
    bool equilibrate_flockwork = true;

    double infection_rate = R0 * recovery_rate / k;

    default_random_engine generator(seed);
    uniform_real_distribution<double> uni_distribution(0.,1.);

    bool start_with_ER = true;


    //start with ER network
    for (size_t i=0; i<N-1; i++)
        for (size_t j=i+1; j<N; j++)
        {
            if (start_with_ER)
            {
                if (uni_distribution(generator)<p)
                    E.push_back( make_pair(i,j) );
            } else {
                E.push_back( make_pair(i,j) );
            }
        }


    SIR_result result;

    result = SIRS(E,N,Q,
                 t_run_total,
                 infection_rate,
                 recovery_rate,
                 susceptible_rate,
                 rewiring_rate,
                 number_of_vaccinated,
                 number_of_infected,
                 use_random_rewiring,
                 equilibrate_flockwork,
                 seed
             );


    for(auto p: result.I_of_t)
    {
        cout << p.first << " " << p.second << endl;
    }

    for(auto p: result.R_of_t)
    {
        cout << p.first << " " << p.second << endl;
    }

    return 0;
}
