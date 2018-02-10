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
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "EqFlockwork.h"
#include "SIR.h"
#include "SIRS.h"
#include "SIS.h"
#include "SI.h"
#include "SI_varying.h"
#include "SIS_varying.h"
#include "SIR_varying.h"
#include "SIRS_varying.h"
#include "FW_P_varying.h"
#include "ResultClasses.h"
#include "test_varying_rate.h"
#include "Barrat_model.h"

using namespace std;
namespace py = pybind11;

PYBIND11_PLUGIN(cFlockwork) {
    py::module m("cFlockwork", "Module to perform fast flockwork simulations");
    
    m.def("SI", &SI, "Simulate an SI process on a flockwork given an initial state as an edge list and varying rewiring rates. Returns time and number of infected as well as time and current R0.",
            py::arg("E"),
            py::arg("N"),
            py::arg("Q"),
            py::arg("infection_rate"),
            py::arg("rewiring_rate") = 1.,
            py::arg("number_of_vaccinated") = 0,
            py::arg("number_of_infected") = 1,
            py::arg("use_random_rewiring") = false,
            py::arg("equilibrate_flockwork") = false,
            py::arg("seed") = 0
            );

    m.def("SI_varying_rates", &SI_varying_rates, "Simulate an SI process on a flockwork given an initial state as an edge list. Returns time and number of infected as well as time and current R0.",
            py::arg("E"),
            py::arg("N"),
            py::arg("Q"),
            py::arg("infection_rate"),
            py::arg("rewiring_rate"),
            py::arg("tmax"),
            py::arg("number_of_vaccinated") = 0,
            py::arg("number_of_infected") = 1,
            py::arg("use_random_rewiring") = false,
            py::arg("equilibrate_flockwork") = false,
            py::arg("seed") = 0
            );

    m.def("SIS", &SIS, "Simulate an SIS process on a flockwork given an initial state as an edge list. Returns time and number of infected as well as time and current R0.",
            py::arg("E"),
            py::arg("N"),
            py::arg("Q"),
            py::arg("t_run_total"),
            py::arg("infection_rate"),
            py::arg("recovery_rate"),
            py::arg("rewiring_rate") = 1.,
            py::arg("number_of_vaccinated") = 0,
            py::arg("number_of_infected") = 1,
            py::arg("use_random_rewiring") = false,
            py::arg("equilibrate_flockwork") = false,
            py::arg("seed") = 0
            );

    m.def("SIS_P", &SIS_P, "Simulate an SIS process on a flockwork P-model given an initial state as an edge list. Returns time and number of infected as well as time and current R0.",
            py::arg("E"),
            py::arg("N"),
            py::arg("P"),
            py::arg("t_run_total"),
            py::arg("infection_rate"),
            py::arg("recovery_rate"),
            py::arg("rewiring_rate") = 1.,
            py::arg("number_of_vaccinated") = 0,
            py::arg("number_of_infected") = 1,
            py::arg("use_random_rewiring") = false,
            py::arg("equilibrate_flockwork") = false,
            py::arg("seed") = 0
            );

    m.def("SIS_varying_rates", &SIS_varying_rates, "Simulate an SIS process on a flockwork given an initial state as an edge list with varying rewiring rate and Q. Returns time and number of infected as well as time and current R0.",
            py::arg("E"),
            py::arg("N"),
            py::arg("Q"),
            py::arg("t_run_total"),
            py::arg("infection_rate"),
            py::arg("recovery_rate"),
            py::arg("rewiring_rate"),
            py::arg("tmax"),
            py::arg("number_of_vaccinated") = 0,
            py::arg("number_of_infected") = 1,
            py::arg("use_random_rewiring") = false,
            py::arg("equilibrate_flockwork") = false,
            py::arg("seed") = 0
            );

    m.def("SIS_varying_rates_P", &SIS_varying_rates_P, "Simulate an SIS process on a flockwork P-model given an initial state as an edge list with varying rewiring rate and P. Returns time and number of infected as well as time and current R0.",
            py::arg("E"),
            py::arg("N"),
            py::arg("P"),
            py::arg("t_run_total"),
            py::arg("infection_rate"),
            py::arg("recovery_rate"),
            py::arg("rewiring_rate"),
            py::arg("tmax"),
            py::arg("number_of_vaccinated") = 0,
            py::arg("number_of_infected") = 1,
            py::arg("use_random_rewiring") = false,
            py::arg("equilibrate_flockwork") = false,
            py::arg("seed") = 0
            );

    m.def("SIS_varying_rates_P_neighbor_affinity", &SIS_varying_rates_P_neighbor_affinity, "Simulate an SIS process on a flockwork P-model given an initial state as an edge list with varying rewiring rate and P. Accounts for the integrated social structure of the network. Returns time and number of infected as well as time and current R0.",
            py::arg("E"),
            py::arg("N"),
            py::arg("P"),
            py::arg("t_run_total"),
            py::arg("infection_rate"),
            py::arg("recovery_rate"),
            py::arg("rewiring_rate"),
            py::arg("neighbor_affinity"),
            py::arg("tmax"),
            py::arg("number_of_vaccinated") = 0,
            py::arg("number_of_infected") = 1,
            py::arg("use_random_rewiring") = false,
            py::arg("equilibrate_flockwork") = false,
            py::arg("use_preferential_node_selection") = false,
            py::arg("use_unweighted_k_for_selection") = false,
            py::arg("seed") = 0
            );

    m.def("SIS_varying_rates_P_neighbor_affinity_infected_nodes", &SIS_varying_rates_P_neighbor_affinity_infected_nodes, "Simulate an SIS process on a flockwork P-model with varying rewiring rate and P. Accounts for the accum. social network. Returns time and number of infected as well as time and current R0.",
            py::arg("E"),
            py::arg("N"),
            py::arg("P"),
            py::arg("t_run_total"),
            py::arg("infection_rate"),
            py::arg("recovery_rate"),
            py::arg("rewiring_rate"),
            py::arg("neighbor_affinity"),
            py::arg("tmax"),
            py::arg("vaccinated_nodes"),
            py::arg("infected_nodes"),
            py::arg("use_random_rewiring") = false,
            py::arg("equilibrate_flockwork") = false,
            py::arg("use_preferential_node_selection") = false,
            py::arg("use_unweighted_k_for_selection") = false,
            py::arg("seed") = 0
            );

    m.def("SIR", &SIR, "Simulate an SIR process on a flockwork given an initial state as an edge list. Returns time and number of infected as well as time and current R0.",
            py::arg("E"),
            py::arg("N"),
            py::arg("Q"),
            py::arg("infection_rate"),
            py::arg("recovery_rate"),
            py::arg("rewiring_rate") = 1.,
            py::arg("t_run_total") = 0,
            py::arg("number_of_vaccinated") = 0,
            py::arg("number_of_infected") = 1,
            py::arg("use_random_rewiring") = false,
            py::arg("equilibrate_flockwork") = false,
            py::arg("seed") = 0
            );

    m.def("SIR_varying_rates", &SIR_varying_rates, "Simulate an SIR process on a flockwork given an initial state as an edge list with varying rewiring rate and Q. Returns time and number of infected as well as time and current R0.",
            py::arg("E"),
            py::arg("N"),
            py::arg("Q"),
            py::arg("infection_rate"),
            py::arg("recovery_rate"),
            py::arg("rewiring_rate"),
            py::arg("tmax"),
            py::arg("t_run_total") = 0,
            py::arg("number_of_vaccinated") = 0,
            py::arg("number_of_infected") = 1,
            py::arg("use_random_rewiring") = false,
            py::arg("equilibrate_flockwork") = false,
            py::arg("seed") = 0
            );

    m.def("flockwork_P_varying_rates", &flockwork_P_varying_rates, "Simulate a flockwork P-model given an initial state as an edge list with varying rewiring rate and varying P. Returns time points and concurrent edge changes.",
            py::arg("E"),
            py::arg("N"),
            py::arg("P"),
            py::arg("t_run_total"),
            py::arg("rewiring_rate"),
            py::arg("tmax"),
            py::arg("use_random_rewiring") = false,
            py::arg("equilibrate_flockwork") = false,
            py::arg("seed") = 0
         );

    m.def("ZSBB_model", &ZSBB_model, "Simulate model after Zhao, Stehle, Bianconi, and Barrat.",
            py::arg("E"),
            py::arg("N"),
            py::arg("lambda"),
            py::arg("b0"),
            py::arg("b1"),
            py::arg("t_run_total"),
            py::arg("seed") = 0,
            py::arg("verbose") = false
         );

    m.def("flockwork_P_varying_rates_neighbor_affinity", &flockwork_P_varying_rates_neighbor_affinity, "Simulate a flockwork P-model given an initial state as an edge list with varying rewiring rate and varying P. Rewiring neighbors are chosen according to a neighbor affinity value. Returns time points and concurrent edge changes.",
            py::arg("E"),
            py::arg("N"),
            py::arg("P"),
            py::arg("t_run_total"),
            py::arg("rewiring_rate"),
            py::arg("neighbor_affinity"),
            py::arg("tmax"),
            py::arg("use_random_rewiring") = false,
            py::arg("equilibrate_flockwork") = false,
            py::arg("use_preferential_node_selection") = false,
            py::arg("use_unweighted_k_for_selection") = false,
            py::arg("seed") = 0
         );

    m.def("SIRS", &SIRS, "Simulate an SIRS process on a flockwork given an initial state as an edge list. Returns time and number of infected as well as time and current R0.",
            py::arg("E"),
            py::arg("N"),
            py::arg("Q"),
            py::arg("t_run_total"),
            py::arg("infection_rate"),
            py::arg("recovery_rate"),
            py::arg("susceptible_rate"),
            py::arg("rewiring_rate") = 1.,
            py::arg("number_of_vaccinated") = 0,
            py::arg("number_of_infected") = 1,
            py::arg("use_random_rewiring") = false,
            py::arg("equilibrate_flockwork") = false,
            py::arg("seed") = 0
            );

    m.def("SIRS_varying_rates", &SIRS_varying_rates, "Simulate an SIRS process on a flockwork given an initial state as an edge list with varying rewiring rates and Q. Returns time and number of infected as well as time and current R0.",
            py::arg("E"),
            py::arg("N"),
            py::arg("Q"),
            py::arg("t_run_total"),
            py::arg("infection_rate"),
            py::arg("recovery_rate"),
            py::arg("susceptible_rate"),
            py::arg("rewiring_rate"),
            py::arg("tmax"),
            py::arg("number_of_vaccinated") = 0,
            py::arg("number_of_infected") = 1,
            py::arg("use_random_rewiring") = false,
            py::arg("equilibrate_flockwork") = false,
            py::arg("seed") = 0
            );

    m.def("equilibrate", &equilibrate_edgelist_seed, "Equilibrates a flockwork given an initial state as an edge list. Returns new edge list.",
            py::arg("E"),
            py::arg("N"),
            py::arg("Q"),
            py::arg("seed"),
            py::arg("t_max") = 0,
            py::arg("use_Q_as_P") = false
            );

    m.def("equilibrate_P", &equilibrate_edgelist_seed, "Equilibrates a flockwork P-model given an initial state as an edge list. Returns new edge list.",
            py::arg("E"),
            py::arg("N"),
            py::arg("P"),
            py::arg("seed"),
            py::arg("t_max") = 0,
            py::arg("use_Q_as_P") = true
            );

    m.def("simulate", &simulate_flockwork, "Simulates a flockwork given an initial state as an edge list. Returns new edge list.",
            py::arg("E"),
            py::arg("N"),
            py::arg("Q"),
            py::arg("seed"),
            py::arg("num_timesteps")
            );

    m.def("simulate_P", &simulate_flockwork_P, "Simulates a flockwork P-model given an initial state as an edge list. Returns new edge list.",
            py::arg("E"),
            py::arg("N"),
            py::arg("P"),
            py::arg("seed"),
            py::arg("num_timesteps")
            );

    m.def("simulate_P_group_lifetimes", &simulate_flockwork_P_group_life_time, "Simulates a flockwork P-model with initial state of all nodes unconnected. Returns a list of group life times.",
            py::arg("N"),
            py::arg("P"),
            py::arg("seed"),
            py::arg("num_timesteps")
            );

    m.def("gillespie_tau_and_event_varying_gamma",&gillespie_tau_and_event_varying_gamma,"Get a Gillespie tau and event number for varying rewiring rates",
            py::arg("standard_rates"),
            py::arg("gamma"),
            py::arg("t0"),
            py::arg("ti"),
            py::arg("t_max"),
            py::arg("seed")
         );


    py::class_<SIR_result>(m,"SIR_result")
        .def(py::init<>())
        .def_readwrite("I_of_t", &SIR_result::I_of_t)
        .def_readwrite("R_of_t", &SIR_result::R_of_t)
        .def_readwrite("SI_of_t", &SIR_result::SI_of_t)
        .def_readwrite("R0_of_t", &SIR_result::R0_of_t)
        .def_readwrite("edge_list", &SIR_result::edge_list);

    py::class_<SIS_result>(m,"SIS_result")
        .def(py::init<>())
        .def_readwrite("I_of_t", &SIS_result::I_of_t)
        .def_readwrite("SI_of_t", &SIS_result::SI_of_t)
        .def_readwrite("R0_of_t", &SIS_result::R0_of_t)
        .def_readwrite("edge_list", &SIS_result::edge_list);

    py::class_<edge_changes>(m,"edge_changes")
        .def(py::init<>())
        .def_readwrite("t", &edge_changes::t)
        .def_readwrite("edges_out", &edge_changes::edges_out)
        .def_readwrite("edges_in", &edge_changes::edges_in);


    return m.ptr();

}
