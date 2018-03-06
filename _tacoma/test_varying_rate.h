#ifndef __TEST_VARYING_H__
#define __TEST_VARYING_H__

#include "Utilities.h"
#include <random>
#include <tuple>

using namespace std;

tuple < int, double, int > gillespie_tau_and_event_varying_gamma( 
                vector < double > & standard_rates,
                vector < pair < double, double > > & gamma,
                double t0,
                size_t ti,
                double t_max,
                size_t seed
             );

#endif
