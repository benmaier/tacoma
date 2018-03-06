#include "Utilities.h"
#include "test_varying_rate.h"
#include <random>
#include <ctime>
#include <tuple>


using namespace std;

tuple < int, double, int > gillespie_tau_and_event_varying_gamma( 
                vector < double > & standard_rates,
                vector < pair < double, double > > & gamma,
                double t0,
                size_t ti,
                double t_max,
                size_t seed
               )
{
    //initialize random generators
    default_random_engine generator(seed);
    uniform_real_distribution<double> uni_distribution(0.,1.);

    double tau = 0.0;
    size_t event = 0;

    get_gillespie_tau_and_event_with_varying_gamma(
                standard_rates,
                gamma,
                t0,
                t_max,
                ti,
                tau,
                event,
                generator, 
                uni_distribution
                );

    return make_tuple(ti,tau,event);
}
