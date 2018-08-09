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
#include "Flockwork.h"
#include "FW_P_varying.h"
#include "FW_P_varying_alpha_beta.h"
#include "ResultClasses.h"
#include "test_varying_rate.h"
#include "ZSBB_model.h"
#include "dyn_RGG.h"
#include "measurements.h"
#include "resampling.h"
#include "social_trajectories.h"
#include "verify_formats.h"
#include "SIS.h"
#include "SIR.h"
#include "SI.h"
#include "SIRS.h"
#include "dyn_gillespie.h"
#include "conversion.h"
#include "concatenation.h"
#include "flockwork_parameter_estimation.h"

using namespace std;
namespace py = pybind11;

PYBIND11_MODULE(_tacoma, m) {
    //TemporAl COntact Modeling and Analysis. Provides fast tools to analyze temporal contact networks and simulate Gillespie processes on them.
    m.doc() = R"pbdoc(
        TemporAl COntact Modeling and Analysis. Provides fast tools to analyze temporal contact networks and simulate Gillespie processes on them.

        .. currentmodule:: _tacoma

        Temporal network classes
        ------------------------------

        .. autosummary::
            :toctree: _generate

            edge_lists
            edge_changes
            edge_trajectories

        Compartmental infection models
        ------------------------------

        .. autosummary::
            :toctree: _generate

            SI
            SIS
            SIR
            SIRS

        Analysis classes
        ---------------

        .. autosummary::
            :toctree: _generate

            SI

        Helper classes
        --------------

        .. autosummary::
            :toctree: _generate

            edge_trajectory_entry
            social_trajectory_entry
            edge_weight
    )pbdoc";
    //m.attr("__name__") = "tacoma.core";
    
    m.def("flockwork_P_varying_rates", &flockwork_P_varying_rates, "Simulate a flockwork P-model given an initial state as an edge list with varying rewiring rate and varying P. Returns time points and concurrent edge changes.",
            py::arg("E"),
            py::arg("N"),
            py::arg("P"),
            py::arg("t_run_total"),
            py::arg("rewiring_rate"),
            py::arg("tmax"),
            py::arg("use_random_rewiring") = false,
            py::arg("seed") = 0
         );

    m.def("flockwork_P_varying_rates_for_each_node", &flockwork_P_varying_rates_for_each_node, 
            "Simulate a flockwork P-model given an initial state as an edge list with varying rewiring rates and varying P (varying both over time and for each node). Returns time points and concurrent edge changes.",
            py::arg("E"),
            py::arg("N"),
            py::arg("P"),
            py::arg("t_run_total"),
            py::arg("rewiring_rates"),
            py::arg("tmax"),
            py::arg("use_random_rewiring") = false,
            py::arg("seed") = 0
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
            py::arg("use_preferential_node_selection") = false,
            py::arg("use_unweighted_k_for_selection") = false,
            py::arg("seed") = 0
         );

    m.def("flockwork_alpha_beta_varying_rates", &flockwork_alpha_beta_varying_rates, "Simulate a flockwork alpha-beta-model given an initial state as an edge list with varying reconnection rate alpha and varying disconnection rate beta. Returns time points and concurrent edge changes.",
            py::arg("E"),
            py::arg("N"),
            py::arg("reconnection_rate"),
            py::arg("disconnection_rate"),
            py::arg("t_run_total"),
            py::arg("tmax"),
            py::arg("use_random_rewiring") = false,
            py::arg("seed") = 0
         );

    m.def("flockwork_P_varying_rates_for_each_node", &flockwork_alpha_beta_varying_rates_for_each_node, 
            "Simulate a flockwork alpha-beta-model given an initial state as an edge list with varying reconnection rate alpha and varying disconnection rate beta (varying both over time and for each node). Returns time points and concurrent edge changes.",
            py::arg("E"),
            py::arg("N"),
            py::arg("reconnection_rates"),
            py::arg("disconnection_rates"),
            py::arg("t_run_total"),
            py::arg("tmax"),
            py::arg("use_random_rewiring") = false,
            py::arg("seed") = 0
         );

    m.def("equilibrate_flockwork_Q", &equilibrate_edgelist_seed, "Equilibrates a flockwork given an initial state as an edge list. Returns new edge list.",
            py::arg("E"),
            py::arg("N"),
            py::arg("Q"),
            py::arg("seed"),
            py::arg("t_max") = 0,
            py::arg("use_Q_as_P") = false
            );

    m.def("equilibrate_flockwork_P", &equilibrate_edgelist_seed, "Equilibrates a flockwork P-model given an initial state as an edge list. Returns new edge list.",
            py::arg("E"),
            py::arg("N"),
            py::arg("P"),
            py::arg("seed"),
            py::arg("t_max") = 0,
            py::arg("use_Q_as_P") = true
            );

    m.def("simulate_flockwork_Q", &simulate_flockwork, "Simulates a flockwork given an initial state as an edge list. Returns new edge list.",
            py::arg("E"),
            py::arg("N"),
            py::arg("Q"),
            py::arg("seed"),
            py::arg("num_timesteps")
            );

    m.def("simulate_flockwork_P", &simulate_flockwork_P, "Simulates a flockwork P-model given an initial state as an edge list. Returns new edge list.",
            py::arg("E"),
            py::arg("N"),
            py::arg("P"),
            py::arg("seed"),
            py::arg("num_timesteps")
            );

    m.def("simulate_flockwork_P_group_lifetimes", &simulate_flockwork_P_group_life_time, "Simulates a flockwork P-model with initial state of all nodes unconnected. Returns a list of group life times.",
            py::arg("N"),
            py::arg("P"),
            py::arg("seed"),
            py::arg("num_timesteps")
            );

    m.def("measure_group_sizes_and_durations_for_edge_lists", &measure_group_sizes_and_durations, "Get a temporal network as a list of edge lists, list of times and number of nodes N and return a list of contact durations, a list of group size histograms and a list of durations lists, one for each group size.",
            py::arg("edge_lists"),
            py::arg("ignore_size_histograms") = false,
            py::arg("verbose") = false
         );

    m.def("measure_group_sizes_and_durations_for_edge_changes", &measure_group_sizes_and_durations_for_edge_changes, "Get a temporal network as a list of edge lists, list of times and number of nodes N and return a list of contact durations, a list of group size histograms and a list of durations lists, one for each group size.",
            py::arg("edge_changes"),
            py::arg("ignore_size_histogram_differences") = false,
            py::arg("verbose") = false
         );

    m.def("mean_degree_from_edge_lists", &mean_degree_from_edge_lists, "Get a list of pairs, where the first entry is a time point and the second entry is the mean degree of the current state, given an `edge_lists` instance.",
            py::arg("edge_lists")
         );

    m.def("mean_degree_from_edge_changes", &mean_degree_from_edge_changes, "Get a list of pairs, where the first entry is a time point and the second entry is the mean degree of the current state, given an `edge_changes` instance.",
            py::arg("edge_changes")
         );

    m.def("degree_distribution_from_edge_lists", &degree_distribution_from_edge_lists, "Get a list of doubles, where the k-th entry of the list is the time-averaged probability that a node has degree k.",
            py::arg("edge_lists")
         );

    m.def("degree_distribution_from_edge_changes", &degree_distribution_from_edge_changes, "Get a list of doubles, where the k-th entry of the list is the time-averaged probability that a node has degree k.",
            py::arg("edge_changes")
         );

    m.def("get_edge_counts", &get_edge_counts, "Given an instance of `edge_changes`, returns the lists `m_in`, `m_out` and `m`, which count the edge events as well as the edges at the given times. Note that `m` contains one entry more than `m_in` and `m_out` due to the initial edge list.",
            py::arg("edge_changes")
         );

    m.def("sample_from_edge_lists", &sample_from_edge_lists, "Get a temporal network as a list of edge lists, given another instance of `edge_lists`, but resampled every dt. Alternatively, provide a number of time steps to divide (tmax-t0) into. if `sample_aggregates` is `True`, this does not sample the network state after each dt, but rather gives a graph of all edges being present in the last time step.",
            py::arg("edge_lists"),
            py::arg("dt") = 0.0,
            py::arg("N_time_steps") = 0,
            py::arg("sample_aggregates") = false,
            py::arg("verbose") = false
         );

    m.def("sample_from_edge_changes", &sample_from_edge_changes, "Get a temporal network as a list of edge lists, given an instance of `edge_changes`, but resampled every dt. Alternatively, provide a number of time steps to divide (tmax-t0) into. if `sample_aggregates` is `True`, this does not sample the network state after each dt, but rather gives a graph of all edges being present in the last time step.",
            py::arg("edge_changes"),
            py::arg("dt") = 0.0,
            py::arg("N_time_steps") = 0,
            py::arg("sample_aggregates") = false,
            py::arg("verbose") = false
         );

    m.def("bin_from_edge_lists", &bin_from_edge_lists, "Get a temporal network as a list of edge lists, given another instance of `edge_lists`, but binned for every dt. Alternatively, provide a number of time steps to divide (tmax-t0) into. This does not sample the network state after each dt, but rather gives a graph of all edges being present in the last time step.",
            py::arg("edge_lists"),
            py::arg("dt") = 0.0,
            py::arg("N_time_steps") = 0,
            py::arg("verbose") = false
         );

    m.def("bin_from_edge_changes", &bin_from_edge_changes, "Get a temporal network as a list of edge lists, given an instance of `edge_changes`, but binned every dt. Alternatively, provide a number of time steps to divide (tmax-t0) into. This does not sample the network state after each dt, but rather gives a graph of all edges being present in the last time step.",
            py::arg("edge_changes"),
            py::arg("dt") = 0.0,
            py::arg("N_time_steps") = 0,
            py::arg("verbose") = false
         );

    m.def("binned_social_trajectory_from_edge_changes", &binned_social_trajectory_from_edge_changes, "Get a social trajectory of a node, given an instance of `edge_changes`, binned every dt (list of group indices this node was part of in [t,t+dt) ).",
            py::arg("edge_changes"),
            py::arg("node"),
            py::arg("dt") = 0.0,
            py::arg("N_time_steps") = 0,
            py::arg("verbose") = false
         );

    m.def("social_trajectory_from_edge_changes", &social_trajectory_from_edge_changes, "Get a social trajectory of a node, given an instance of `edge_changes`.",
            py::arg("edge_changes"),
            py::arg("node"),
            py::arg("verbose") = false
         );

    m.def("social_trajectory_from_edge_lists", &social_trajectory_from_edge_lists, "Get a social trajectory of a node, given an instance of `edge_lists`.",
            py::arg("edge_lists"),
            py::arg("node"),
            py::arg("verbose") = false
         );

    m.def("edge_trajectories_from_edge_lists", &edge_trajectories_from_edge_lists, "For each edge, get all time pairs where the edge existed, given an instance of `edge_lists`.",
            py::arg("edge_lists"),
            py::arg("return_edge_similarities") = false,
            py::arg("verbose") = false
         );

    m.def("edge_trajectories_from_edge_changes", &edge_trajectories_from_edge_changes, "For each edge, get all time pairs where the edge existed, given an instance of `edge_changes`.",
            py::arg("edge_changes"),
            py::arg("return_edge_similarities") = false,
            py::arg("verbose") = false
         );

    m.def("verify_edge_lists", &verify_edge_lists, "For each edge, get all time pairs where the edge existed, given an instance of `edge_lists`.",
            py::arg("edge_lists"),
            py::arg("verbose") = false
         );

    m.def("verify_edge_changes", &verify_edge_changes, "For each edge, get all time pairs where the edge existed, given an instance of `edge_changes`.",
            py::arg("edge_changes"),
            py::arg("verbose") = false
         );

    m.def("convert_edge_lists", &convert_edge_lists, "Convert an instance of `edge_lists` to an instance of `edge_changes`.",
            py::arg("edge_lists"),
            py::arg("verbose") = false
         );

    m.def("convert_edge_changes", &convert_edge_changes, "Convert an instance of `edge_changes` to an instance of `edge_lists`.",
            py::arg("edge_changes"),
            py::arg("verbose") = false
         );

    m.def("concatenate_edge_lists", &concatenate_edge_lists, "Concatenate a list of `edge_lists` to a single instance of `edge_lists`.",
            py::arg("list_of_edge_lists"),
            py::arg("verbose") = false
         );

    m.def("concatenate_edge_changes", &concatenate_edge_changes, "Convert a list of `edge_changes` to a single instance of `edge_changes`.",
            py::arg("list_of_edge_changes"),
            py::arg("verbose") = false
         );

    m.def("binned_social_trajectory_from_edge_lists", &binned_social_trajectory_from_edge_lists, "Get a social trajectory of a node, given an instance of `edge_lists`, binned every dt (list of group indices this node was part of in [t,t+dt) ).",
            py::arg("edge_lists"),
            py::arg("node"),
            py::arg("dt") = 0.0,
            py::arg("N_time_steps") = 0,
            py::arg("verbose") = false
         );

    m.def("get_flockwork_P_args", &get_flockwork_P_args, "Calculate the rewiring_rate gamma(t) and probability to stay alone P(t) as well as the other important parameters to simulate a flockwork_P model with varying rates.",
            py::arg("edge_changes"),
            py::arg("dt") = 0.0,
            py::arg("N_time_steps") = 0,
            py::arg("k_over_k_real_scaling") = 1.0,
            py::arg("gamma_scaling") = 1.0,
            py::arg("P_scaling") = 1.0,
            py::arg("aggregated_network") = map < pair < size_t, size_t >, double >(),
            py::arg("ensure_empty_network") = false,
            py::arg("adjust_last_bin_if_dt_does_not_fit") = false,
            py::arg("verbose") = false
         );

    m.def("get_flockwork_P_node_parameters_gamma_and_P", get_node_gamma_and_P,
            R"pbdoc(Calculate the mean node specific activity gamma_i and reconnection probability P_i for each node.)pbdoc",
            py::arg("edge_changes"),
            py::arg("gamma"),
            py::arg("P"),
            py::arg("use_event_rate_method") = false
        );

    m.def("get_flockwork_P_node_parameters_alpha_and_beta", get_node_gamma_and_P,
            R"pbdoc(Calculate the mean node specific activity gamma_i and reconnection probability P_i for each node.)pbdoc",
            py::arg("edge_changes"),
            py::arg("gamma"),
            py::arg("P"),
            py::arg("use_event_rate_method") = true
        );

    m.def("estimate_k_scaling_gradient_descent", &estimate_k_scaling_gradient_descent, "Estimate the scaling of <k> that's necessary to revert the effects of binning using gradient descent",
            py::arg("edge_changes"),
            py::arg("dt_for_inference"),
            py::arg("dt_for_binning"),
            py::arg("measurements_per_configuration") = 6,
            py::arg("learning_rate") = 0.5,
            py::arg("relative_error") = 1e-2,
            py::arg("N_eval_max") = 100,
            py::arg("verbose") = true
        );

    m.def("estimate_k_scaling_gradient_descent_RMSE", &estimate_k_scaling_gradient_descent_RMSE, "Estimate the scaling of <k> that's necessary to revert the effects of binning using gradient descent",
            py::arg("edge_changes"),
            py::arg("dt_for_inference"),
            py::arg("dt_for_binning"),
            py::arg("measurements_per_configuration") = 6,
            py::arg("learning_rate") = 0.5,
            py::arg("relative_error") = 1e-2,
            py::arg("N_eval_max") = 100,
            py::arg("verbose") = true
        );

    m.def("ZSBB_model", &ZSBB_model, "Simulate model after Zhao, Stehle, Bianconi, and Barrat.",
            py::arg("E"),
            py::arg("N"),
            py::arg("lambda"),
            py::arg("b0"),
            py::arg("b1"),
            py::arg("t_run_total") = 0,
            py::arg("max_edge_events_to_end_simulation") = 0,
            py::arg("t_equilibration") = 0,
            py::arg("seed") = 0,
            py::arg("record_sizes_and_durations") = false,
            py::arg("return_after_equilibration_only") = false,
            py::arg("verbose") = false
         );

    m.def("dynamic_RGG", &dynamic_RGG, "Simulate dynamic random geometric graph model.",
            py::arg("N"),
            py::arg("t_run_total"),
            py::arg("step_distance") = 0.0,
            py::arg("mean_link_duration") = 0.0,
            py::arg("critical_density") = 0.65,
            py::arg("periodic_boundary_conditions_for_link_building") = false,
            py::arg("record_sizes_and_durations") = false,
            py::arg("seed") = 0,
            py::arg("verbose") = false
         );

    m.def("gillespie_SIS_on_edge_lists",&gillespie_on_edge_lists<SIS>,"Perform a Gillespie SIS simulation on edge lists. Needs an instance of tacoma.edge_lists and an instance of tacoma.SIS.",
            py::arg("edge_lists"),
            py::arg("SIS"),
            py::arg("is_static") = false,
            py::arg("verbose") = false
            );

    m.def("gillespie_SIS_on_edge_changes",&gillespie_on_edge_changes<SIS>,"Perform a Gillespie SIS simulation on edge changes. Needs an instance of tacoma.edge_changes and an instance of tacoma.SIS.",
            py::arg("edge_changes"),
            py::arg("SIS"),
            py::arg("verbose") = false
            );

    m.def("gillespie_SI_on_edge_lists",&gillespie_on_edge_lists<SI>,"Perform a Gillespie SI simulation on edge lists. Needs an instance of tacoma.edge_lists and an instance of tacoma.SI.",
            py::arg("edge_lists"),
            py::arg("SI"),
            py::arg("is_static") = false,
            py::arg("verbose") = false
            );

    m.def("gillespie_SI_on_edge_changes",&gillespie_on_edge_changes<SI>,"Perform a Gillespie SI simulation on edge changes. Needs an instance of tacoma.edge_changes and an instance of tacoma.SI.",
            py::arg("edge_changes"),
            py::arg("SI"),
            py::arg("verbose") = false
            );

    m.def("gillespie_SIR_on_edge_lists",&gillespie_on_edge_lists<SIR>,"Perform a Gillespie SIR simulation on edge lists. Needs an instance of tacoma.edge_lists and an instance of tacoma.SIR.",
            py::arg("edge_lists"),
            py::arg("SIR"),
            py::arg("is_static") = false,
            py::arg("verbose") = false
            );

    m.def("gillespie_SIR_on_edge_changes",&gillespie_on_edge_changes<SIR>,"Perform a Gillespie SIR simulation on edge changes. Needs an instance of tacoma.edge_changes and an instance of tacoma.SIR.",
            py::arg("edge_changes"),
            py::arg("SIR"),
            py::arg("verbose") = false
            );

    m.def("gillespie_SIRS_on_edge_lists",&gillespie_on_edge_lists<SIRS>,"Perform a Gillespie SIRS simulation on edge lists. Needs an instance of tacoma.edge_lists and an instance of tacoma.SIRS.",
            py::arg("edge_lists"),
            py::arg("SIRS"),
            py::arg("is_static") = false,
            py::arg("verbose") = false
            );

    m.def("gillespie_SIRS_on_edge_changes",&gillespie_on_edge_changes<SIRS>,"Perform a Gillespie SIRS simulation on edge changes. Needs an instance of tacoma.edge_changes and an instance of tacoma.SIRS.",
            py::arg("edge_changes"),
            py::arg("SIRS"),
            py::arg("verbose") = false
            );

    py::class_<edge_changes>(m,"edge_changes", R"pbdoc(Description of a temporal network by listing the changes of edges at a certain time.)pbdoc")
        .def(py::init<>())
        .def(py::init<const edge_changes_with_histograms &>(),
                py::arg("edge_changes_with_histograms"),
                "Initialize from an instance of :mod:`edge_changes_with_histograms`"
            )
        .def("copy_from", &edge_changes::copy_from, R"pbdoc([deprecated] copy the relevant information from an instance of :mod:`edge_changes_with_histograms`)pbdoc")
        .def_readwrite("int_to_node", &edge_changes::int_to_node, R"pbdoc(A dictionary int -> string which keeps the original node names.)pbdoc")
        .def_readwrite("time_unit", &edge_changes::time_unit, R"pbdoc(A string containing the unit of time for this network.)pbdoc")
        .def_readwrite("notes", &edge_changes::notes, R"pbdoc(A string containing additional notes for this network.)pbdoc")
        .def_readwrite("t", &edge_changes::t, R"pbdoc(An ordered list containing the time points at which changes occur.)pbdoc")
        .def_readwrite("edges_out", &edge_changes::edges_out, R"pbdoc(A list containing the edges leaving the network at the correspoding time in `t`)pbdoc")
        .def_readwrite("edges_in", &edge_changes::edges_in, R"pbdoc(A list containing the edges coming in to the network at the correspoding time in `t`)pbdoc")
        .def_readwrite("N", &edge_changes::N, R"pbdoc(Number of nodes)pbdoc")
        .def_readwrite("t0", &edge_changes::t0, R"pbdoc(The initial time)pbdoc")
        .def_readwrite("tmax", &edge_changes::tmax, R"pbdoc(The final time)pbdoc")
        .def_readwrite("edges_initial", &edge_changes::edges_initial, R"pbdoc(A list containing the edges of the network at time `t0`)pbdoc");

    py::class_<edge_lists>(m,"edge_lists", R"pbdoc(
            Description of a temporal network by listing 
            the edges of the networks at a certain time.
            )pbdoc")
        .def(py::init<>())
        .def(py::init<const edge_lists_with_histograms &>(),
                py::arg("edge_lists_with_histograms"),
                "Initialize from an instance of :mod:`edge_lists_with_histograms`"
            )
        .def("copy_from", &edge_lists::copy_from, R"pbdoc([deprecated] copy the relevant information from an instance of :mod:`edge_lists_with_histograms`)pbdoc")
        .def_readwrite("int_to_node", &edge_lists::int_to_node, R"pbdoc(A dictionary int -> string which keeps the original node names.)pbdoc")
        .def_readwrite("time_unit", &edge_lists::time_unit, R"pbdoc(A string containing the unit of time for this network.)pbdoc")
        .def_readwrite("notes", &edge_lists::notes, R"pbdoc(A string containing additional notes for this network.)pbdoc")
        .def_readwrite("t", &edge_lists::t, R"pbdoc(An ordered list containing the time points at which the new edge list becomes active.)pbdoc")
        .def_readwrite("edges", &edge_lists::edges, R"pbdoc(A list containing the edge list of the network at the correspoding time in `t`)pbdoc")
        .def_readwrite("N", &edge_lists::N, R"pbdoc(Number of nodes)pbdoc")
        .def_readwrite("tmax", &edge_lists::tmax, R"pbdoc(The final time)pbdoc");

    py::class_<edge_lists_with_histograms>(m,"edge_lists_with_histograms",R"pbdoc(
            Similar to the :mod:`edge_lists` class but with additional analysis results.
        )pbdoc")
        .def(py::init<>())
        .def_readwrite("t", &edge_lists_with_histograms::t, R"pbdoc(
            An ordered list containing the time points 
            at which the new edge list becomes active.
            )pbdoc")
        .def_readwrite("edges", &edge_lists_with_histograms::edges, R"pbdoc(
                A list containing the edge list of the network at the correspoding time in `t`
        )pbdoc")
        .def_readwrite("size_histograms", &edge_lists_with_histograms::size_histograms)
        .def_readwrite("group_durations", &edge_lists_with_histograms::group_durations)
        .def_readwrite("N", &edge_lists_with_histograms::N)
        .def_readwrite("tmax", &edge_lists_with_histograms::tmax);

    py::class_<edge_changes_with_histograms>(m,"edge_changes_with_histograms")
        .def(py::init<>())
        .def_readwrite("t", &edge_changes_with_histograms::t)
        .def_readwrite("edges_in", &edge_changes_with_histograms::edges_in)
        .def_readwrite("edges_out", &edge_changes_with_histograms::edges_out)
        .def_readwrite("initial_size_histogram", &edge_changes_with_histograms::initial_size_histogram)
        .def_readwrite("group_changes", &edge_changes_with_histograms::group_changes)
        .def_readwrite("final_size_histogram", &edge_changes_with_histograms::final_size_histogram)
        .def_readwrite("contact_durations", &edge_changes_with_histograms::contact_durations)
        .def_readwrite("inter_contact_durations", &edge_changes_with_histograms::inter_contact_durations)
        .def_readwrite("group_durations", &edge_changes_with_histograms::group_durations)
        .def_readwrite("N", &edge_changes_with_histograms::N)
        .def_readwrite("t0", &edge_changes_with_histograms::t0)
        .def_readwrite("tmax", &edge_changes_with_histograms::tmax)
        .def_readwrite("edges_initial", &edge_changes_with_histograms::edges_initial);

    py::class_<group_sizes_and_durations>(m,"group_sizes_and_durations")
        .def(py::init<>())
        .def_readwrite("contact_durations", &group_sizes_and_durations::contact_durations)
        .def_readwrite("size_histograms", &group_sizes_and_durations::size_histograms)
        .def_readwrite("size_histogram_differences", &group_sizes_and_durations::size_histogram_differences)
        .def_readwrite("group_durations", &group_sizes_and_durations::group_durations)
        .def_readwrite("aggregated_size_histogram", &group_sizes_and_durations::aggregated_size_histogram)
        .def_readwrite("aggregated_network", &group_sizes_and_durations::aggregated_network);

    py::class_<edge_weight>(m,"edge_weight",R"pbdoc(Helper class for internal usage. Creates `value = 0` when initiated
    such that one can easily use it as a counter in a map without checking whether or not the object exists.)pbdoc")
        .def(py::init<>())
        .def_readwrite("value", &edge_weight::value);

    py::class_<social_trajectory_entry>(m,"social_trajectory_entry","Each :mod:`social_trajectory_entry` ")
        .def(py::init<>())
        .def_readwrite("hash", &social_trajectory_entry::hash, "A group identifier which has low probability of doubling.")
        .def_readwrite("size", &social_trajectory_entry::size, "Number of nodes within this group.")
        .def_readwrite("time_pairs", &social_trajectory_entry::time_pairs, R"pbdoc(List[Tuple[double, double]] ordered time pairs denoting the time intervals in which the group existed)pbdoc");

    py::class_<edge_trajectory_entry>(m,"edge_trajectory_entry",R"pbdoc(
        This is an entry of an edge-based notation of a temporal 
        network. Instead of getting lists of edges
        ordered in time, a list of :mod:`edge_trajectory_entry` contains, 
        for each edge, a list of time pairs denoting the times the edge exists.
         )pbdoc")
        .def(py::init<>())
        .def_readwrite("edge", &edge_trajectory_entry::edge, R"pbdoc(Tuple[int, int] containing the nodes belonging to this edge.)pbdoc")
        .def_readwrite("time_pairs", &edge_trajectory_entry::time_pairs, R"pbdoc(
            List[Tuple[double, double]] 
            A list containing ordered time pairs :math:`(t_i^{(i)}, t_i^{(f)})`, where each time pair contains 
            the time point :math:`t_i^{(i)}` when the edge is switched on (created) and 
            the time :math:`t_i^{(f)}` when the edge is switched off (deleted).
            )pbdoc"
                );

    py::class_<flockwork_args>(m,"flockwork_args",py::dynamic_attr(), R"pbdoc(
                    The arguments which are passed to the flockwork simulation 
                    function. An instance of this is returned by the parameter estimation procedure.
                )pbdoc")
        .def(py::init<>())
        .def_readwrite("E", &flockwork_args::E, "Edge list of the initial state")
        .def_readwrite("N", &flockwork_args::N, "Number of nodes")
        .def_readwrite("P", &flockwork_args::P, R"pbdoc(A list of floats describing the time-dependent connection probability (corresponding times in `rewiring_rate`))pbdoc")
        .def_readwrite("rewiring_rate", &flockwork_args::rewiring_rate, R"pbdoc(A list of pairs of doubles, each entry contains the time and the rewiring rate per node)pbdoc")
        .def_readwrite("neighbor_affinity", &flockwork_args::neighbor_affinity, R"pbdoc(
        A list, for each node contains two lists, one containing the node's 
        neighbors and the second one containing the node affinity between them.)pbdoc")
        .def_readwrite("tmax", &flockwork_args::tmax, R"pbdoc(The time at which the last value of `P` and the `rewiring_rate` changes (i.e. the maximum time until the parameteres are defined))pbdoc")
        .def_readwrite("m", &flockwork_args::m, "The number of edges in the network at this time")
        .def_readwrite("m_in", &flockwork_args::m_in, "The number of edges being created in the last time interval")
        .def_readwrite("m_out", &flockwork_args::m_out, "The number of edges being deleted in the last time interval")
        .def_readwrite("new_time", &flockwork_args::new_time, "The bin edges of the new time bins")
        ;

    py::class_<edge_trajectories>(m,"edge_trajectories",R"pbdoc(
        Instead of getting lists of edges ordered in time, this description 
        of a temporal network consists of a list of :mod:`edge_trajectory_entry`. Each entry
        contains the edge and a list of time pairs denoting the times the edge exists.

        Optionally, dependent on the function which created this object, this object can
        contain edge similarities, which are defined as follows. 
        Two edges are considered similar when they are connected to the same node at the same time
        Each edge :math:`i`, where
        :math:`i` is the edge's index in the trajectory list, 
    )pbdoc")
        .def(py::init<>())
        .def_readwrite("trajectories", &edge_trajectories::trajectories, R"pbdoc(Each entry of this list has properties `.edge` containing 
            its nodes and `.time_pairs` containing the time intervals when the edge was active.)pbdoc")
        .def_readwrite("edge_similarities", &edge_trajectories::edge_similarities, R"pbdoc(
            Each entry of this list is a triple (i, j, w), where `i` is the `i`-th edge in `trajectories` (similarly for `j`) and `w`
            is their similarity.
        )pbdoc");

    py::class_<SIS>(m,"SIS","Base class for the simulation of an SIS compartmental infection model on a temporal network. Pass this to :mod:`gillespie_SIS` to simulate and retrieve the simulation results.")
        .def(py::init<size_t,double,double,double,size_t,size_t,bool,size_t,bool>(),
                py::arg("N"),
                py::arg("t_simulation"),
                py::arg("infection_rate"),
                py::arg("recovery_rate"),
                py::arg("number_of_initially_infected") = 1, 
                py::arg("number_of_initially_vaccinated") = 0, 
                py::arg("prevent_disease_extinction") = false,
                py::arg("seed") = 0,
                py::arg("verbose") = false,
                R"pbdoc(
                    Parameters
                    ----------
                    N : int
                        Number of nodes in the temporal network.
                    t_simulation : float
                        Maximum time for the simulation to run. Can possibly be greater than the maximum time of the temporal
                        network in which case the temporal network is looped.
                    infection_rate : float
                        Infection rate per :math:`SI`-link (expected number of reaction events :math:`SI\rightarrow II`
                        for a single :math:`SI`-link per dimension of time).
                    recovery_rate : float
                        Recovery rate per infected (expected number of reaction events :math:`I\rightarrow S`
                        for a single infected node per dimension of time).
                    number_of_initially_infected : int, default = 1
                        Number of nodes which will be in the infected compartment at :math:`t = t_0`. Note that the default
                        value 1 is not an ideal initial value as fluctuations may lead to a quick end of the simulation
                        skewing the outcome. I generally recommend to use a number of the order of :math:`N/2`.
                    number_of_initially_vaccinated : int, default = 0
                        Number of nodes which will be in the recovered compartment at :math:`t = t_0`.
                    prevent_disease_extinction : bool, default: False
                        If this is `True`, the recovery of the last infected node will always be prevented.
                    seed : int, default = 0
                        Seed for RNG initialization. If this is 0, the seed will be initialized randomly.
                    verbose : bool, default = False
                        Be talkative.
                )pbdoc"
            )
        .def_readwrite("time", &SIS::time, "A list containing the time points at which one or more of the observables changed.")
    .def_readwrite("R0", &SIS::R0, R"pbdoc(
                   A list containing the basic reproduction number defined as :math:`R_0(t) = \eta\left\langle k \right\rangle(t) / \rho`
                   where :math:`\eta` is the infection rate per link and :math:`\rho` is the recovery rate per node.
                   )pbdoc")
        .def_readwrite("SI", &SIS::SI, "A list containing the number of infected at time :math:`t`.")
        .def_readwrite("I", &SIS::I, "A list containing the number of recovered at time :math:`t`.");

    py::class_<SI>(m,"SI","Base class for the simulation of an SI compartmental infection model on a temporal network. Pass this to :mod:`gillespie_SI` to simulate and retrieve the simulation results.")
        .def(py::init<size_t,double,double,size_t,size_t,size_t,bool>(),
                py::arg("N"),
                py::arg("t_simulation"),
                py::arg("infection_rate"),
                py::arg("number_of_initially_infected") = 1, 
                py::arg("number_of_initially_vaccinated") = 0, 
                py::arg("seed") = 0,
                py::arg("verbose") = false,
                R"pbdoc(
                    Parameters
                    ----------
                    N : int
                        Number of nodes in the temporal network.
                    t_simulation : float
                        Maximum time for the simulation to run. Can possibly be greater than the maximum time of the temporal
                        network in which case the temporal network is looped.
                    infection_rate : float
                        Infection rate per :math:`SI`-link (expected number of reaction events :math:`SI\rightarrow II`
                        for a single :math:`SI`-link per dimension of time).
                    number_of_initially_infected : int, default = 1
                        Number of nodes which will be in the infected compartment at :math:`t = t_0`. 
                    number_of_initially_vaccinated : int, default = 0
                        Number of nodes which will be in the recovered compartment at :math:`t = t_0`.
                    seed : int, default = 0
                        Seed for RNG initialization. If this is 0, the seed will be initialized randomly.
                    verbose : bool, default = False
                        Be talkative.
                )pbdoc"
            )
        .def_readwrite("time", &SI::time, "A list containing the time points at which one or more of the observables changed.")
        .def_readwrite("SI", &SI::_SI, "A list containing the number of infected at time :math:`t`.")
        .def_readwrite("I", &SI::I, "A list containing the number of recovered at time :math:`t`.");

    py::class_<SIR>(m,"SIR","Base class for the simulation of an SIR compartmental infection model on a temporal network. Pass this to :mod:`gillespie_SIR` to simulate and retrieve the simulation results.")
        .def(py::init<size_t,double,double,double,size_t,size_t,size_t,bool>(),
                py::arg("N"),
                py::arg("t_simulation"),
                py::arg("infection_rate"),
                py::arg("recovery_rate"),
                py::arg("number_of_initially_infected") = 1, 
                py::arg("number_of_initially_vaccinated") = 0, 
                py::arg("seed") = 0,
                py::arg("verbose") = false,
                R"pbdoc(
                    Parameters
                    ----------
                    N : int
                        Number of nodes in the temporal network.
                    t_simulation : float
                        Maximum time for the simulation to run. Can possibly be greater than the maximum time of the temporal
                        network in which case the temporal network is looped.
                    infection_rate : float
                        Infection rate per :math:`SI`-link (expected number of reaction events :math:`SI\rightarrow II`
                        for a single :math:`SI`-link per dimension of time).
                    recovery_rate : float
                        Recovery rate per infected (expected number of reaction events :math:`I\rightarrow R`
                        for a single infected node per dimension of time).
                    number_of_initially_infected : int, default = 1
                        Number of nodes which will be in the infected compartment at :math:`t = t_0`. Note that the default
                        value 1 is not an ideal initial value as fluctuations may lead to a quick end of the simulation
                        skewing the outcome. I generally recommend to use a number of the order of :math:`N/2`.
                    number_of_initially_vaccinated : int, default = 0
                        Number of nodes which will be in the recovered compartment at :math:`t = t_0`.
                    seed : int, default = 0
                        Seed for RNG initialization. If this is 0, the seed will be initialized randomly.
                    verbose : bool, default = False
                        Be talkative.
                )pbdoc"
            )
        .def_readwrite("time", &SIR::time, "A list containing the time points at which one or more of the observables changed.")
        .def_readwrite("R0", &SIR::R0, R"pbdoc(
                   A list containing the basic reproduction number defined as :math:`R_0(t) = \eta\left\langle k \right\rangle(t) / \rho`
                   where :math:`\eta` is the infection rate per link and :math:`\rho` is the recovery rate per node.
                   )pbdoc")
        .def_readwrite("SI", &SIR::SI, "A list containing the number of :math:`SI`-links at time :math:`t`.")
        .def_readwrite("I", &SIR::I, "A list containing the number of infected at time :math:`t`.")
        .def_readwrite("R", &SIR::R, "A list containing the number of recovered at time :math:`t`.");

    py::class_<SIRS>(m,"SIRS","Base class for the simulation of an SIRS compartmental infection model on a temporal network. Pass this to :mod:`gillespie_SIRS` to simulate and retrieve the simulation results.")
        .def(py::init<size_t,double,double,double,double,size_t,size_t,size_t,bool>(),
                py::arg("N"),
                py::arg("t_simulation"),
                py::arg("infection_rate"),
                py::arg("recovery_rate"),
                py::arg("waning_immunity_rate"),
                py::arg("number_of_initially_infected") = 1, 
                py::arg("number_of_initially_vaccinated") = 0, 
                py::arg("seed") = 0,
                py::arg("verbose") = false,
                R"pbdoc(
                    Parameters
                    ----------
                    N : int
                        Number of nodes in the temporal network.
                    t_simulation : float
                        Maximum time for the simulation to run. Can possibly be greater than the maximum time of the temporal
                        network in which case the temporal network is looped.
                    infection_rate : float
                        Infection rate per :math:`SI`-link (expected number of reaction events :math:`SI\rightarrow II`
                        for a single :math:`SI`-link per dimension of time).
                    recovery_rate : float
                        Recovery rate per infected (expected number of reaction events :math:`I\rightarrow R`
                        for a single infected node per dimension of time).
                    waning_immunity_rate : float
                        Recovery rate per infected (expected number of reaction events :math:`R\rightarrow S`
                        for a single recovered node per dimension of time).
                    number_of_initially_infected : int, default = 1
                        Number of nodes which will be in the infected compartment at :math:`t = t_0`. Note that the default
                        value 1 is not an ideal initial value as fluctuations may lead to a quick end of the simulation
                        skewing the outcome. I generally recommend to use a number of the order of :math:`N/2`.
                    number_of_initially_vaccinated : int, default = 0
                        Number of nodes which will be in the recovered compartment at :math:`t = t_0`.
                    seed : int, default = 0
                        Seed for RNG initialization. If this is 0, the seed will be initialized randomly.
                    verbose : bool, default = False
                        Be talkative.
                )pbdoc"
            )
        .def_readwrite("time", &SIRS::time, "A list containing the time points at which one or more of the observables changed.")
        .def_readwrite("R0", &SIRS::R0, R"pbdoc(
                    A list containing the basic reproduction number defined as :math:`R_0(t) = \eta\left\langle k \right\rangle(t) / \rho`
                    where :math:`\eta` is the infection rate per link and :math:`\rho` is the recovery rate per node.
            )pbdoc" )
        .def_readwrite("SI", &SIRS::SI, "A list containing the number of :math:`SI`-links at time :math:`t`.")
        .def_readwrite("I", &SIRS::I, "A list containing the number of infected at time :math:`t`.")
        .def_readwrite("R", &SIRS::R, "A list containing the number of recovered at time :math:`t`.");

}
