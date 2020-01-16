/*
 * The MIT License (MIT)
 * Copyright (c) 2016-2018, Benjamin Maier
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
#include "edge_trajectories.h"
#include "verify_formats.h"
#include "SIS.h"
#include "QS_SIS.h"
#include "eSIS.h"
#include "Markov_SIS.h"
#include "markov_dynamics.h"
#include "model_markov.h"
#include "SIS_node_based.h"
#include "coverage_SIS.h"
#include "cluster_size_SIS.h"
#include "SIR.h"
#include "SI.h"
#include "SIRS.h"
#include "dyn_gillespie.h"
#include "model_gillespie.h"
#include "conversion.h"
#include "concatenation.h"
#include "flockwork_parameter_estimation.h"
#include "activity_model.h"
#include "EdgeActivityModel.h"
#include "FlockworkPModel.h"
#include "slice.h"

using namespace std;
namespace py = pybind11;

PYBIND11_MODULE(_tacoma, m)
{
    //TemporAl COntact Modeling and Analysis. Provides fast tools to analyze temporal contact networks and simulate Gillespie processes on them.
    m.doc() = R"pbdoc(
        C++ core of TemporAl COntact Modeling and Analysis. Provides fast tools
        to analyze temporal contact networks and simulate Gillespie processes on them.

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
            node_based_SIS
            eSIS
            QS_SIS
            coverage_SIS
            cluster_size_SIS

        Markov integraion models
        ------------------------

        .. autosummary::
            :toctree: _generate

            MARKOV_SIS

        Temporal network model classes
        ------------------------------

        .. autosummary::
            :toctree: _generate

            EdgeActivityModel
            FlockworkPModel

        Analysis classes
        ----------------

        .. autosummary::
            :toctree: _generate

            edge_lists_with_histograms
            edge_changes_with_histograms
            group_sizes_and_durations

        Helper classes
        --------------

        .. autosummary::
            :toctree: _generate

            edge_trajectory_entry
            social_trajectory_entry
            edge_weight
            EPI
    )pbdoc";

    m.def("flockwork_P_varying_rates",
          &flockwork_P_varying_rates,
          R"pbdoc(Simulate a flockwork P-model given an initial state as an edge list with
          varying rewiring rate and varying reconnection probability P. Returns time points
          and concurrent edge changes.)pbdoc",
          py::arg("E"),
          py::arg("N"),
          py::arg("P"),
          py::arg("t_run_total"),
          py::arg("rewiring_rate"),
          py::arg("tmax"),
          py::arg("use_random_rewiring") = false,
          py::arg("seed") = 0);

    m.def("flockwork_P_varying_rates_for_each_node",
          &flockwork_P_varying_rates_for_each_node,
          R"pbdoc(Simulate a flockwork P-model given an initial state as an edge list with
          varying rewiring rates and varying P (varying both over time and for each node).
          Returns time points and concurrent edge changes.)pbdoc",
          py::arg("E"),
          py::arg("N"),
          py::arg("P"),
          py::arg("t_run_total"),
          py::arg("rewiring_rates"),
          py::arg("tmax"),
          py::arg("use_random_rewiring") = false,
          py::arg("seed") = 0);

    m.def("activity_model",
          &activity_model,
          R"pbdoc(Simulate an activity rate model where active edges become inactive with
          rate omega- and inactive edges become active with rate omega+.)pbdoc",
          py::arg("N"),
          py::arg("rho"),
          py::arg("omega"),
          py::arg("t_run_total"),
          py::arg("seed") = 0,
          py::arg("verbose") = false);

    m.def("activity_model_inefficient",
          &activity_model_inefficient,
          R"pbdoc(Simulate an activity rate model where active edges become inactive with
          rate omega- and inactive edges become active with rate omega+.)pbdoc",
          py::arg("N"),
          py::arg("rho"),
          py::arg("omega"),
          py::arg("t_run_total"),
          py::arg("seed") = 0,
          py::arg("verbose") = false);

    m.def("flockwork_P_varying_rates_neighbor_affinity",
          &flockwork_P_varying_rates_neighbor_affinity,
          R"pbdoc(Simulate a flockwork P-model given an initial state as an edge list with
          varying rewiring rate and varying P. Rewiring neighbors are chosen according to a
          neighbor affinity value. Returns time points and concurrent edge changes.)pbdoc",
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
          py::arg("seed") = 0);

    m.def("flockwork_alpha_beta_varying_rates",
          &flockwork_alpha_beta_varying_rates,
          R"pbdoc(Simulate a flockwork alpha-beta-model given an initial state as an edge
          list with varying reconnection rate alpha and varying disconnection rate beta.
          Returns time points and concurrent edge changes.)pbdoc",
          py::arg("E"),
          py::arg("N"),
          py::arg("reconnection_rate"),
          py::arg("disconnection_rate"),
          py::arg("t_run_total"),
          py::arg("tmax"),
          py::arg("use_random_rewiring") = false,
          py::arg("seed") = 0);

    m.def("flockwork_alpha_beta_varying_rates_with_neighbor_affinity",
          &flockwork_alpha_beta_varying_rates_with_neighbor_affinity,
          R"pbdoc(Simulate a flockwork alpha-beta-model given an initial state as an edge
          list with varying reconnection rate alpha and varying disconnection rate beta.
          Choose neighbors of reacting nodes according to a weighted static network.
          Returns time points and concurrent edge changes.)pbdoc",
          py::arg("E"),
          py::arg("N"),
          py::arg("reconnection_rate"),
          py::arg("disconnection_rate"),
          py::arg("neighbor_affinity"),
          py::arg("t_run_total"),
          py::arg("tmax"),
          py::arg("use_random_rewiring") = false,
          py::arg("seed") = 0);

    m.def("flockwork_alpha_beta_varying_rates_for_each_node",
          &flockwork_alpha_beta_varying_rates_for_each_node,
          R"pbdoc(Simulate a flockwork alpha-beta-model given an initial state as an edge
          list with varying reconnection rate alpha and varying disconnection rate beta
          (varying both over time and for each node). Returns time
          points and concurrent edge changes.)pbdoc",
          py::arg("E"),
          py::arg("N"),
          py::arg("reconnection_rates"),
          py::arg("disconnection_rates"),
          py::arg("t_run_total"),
          py::arg("tmax"),
          py::arg("use_random_rewiring") = false,
          py::arg("seed") = 0);

    m.def("flockwork_alpha_beta_varying_rates_for_each_node_with_neighbor_affinity",
          &flockwork_alpha_beta_varying_rates_for_each_node_with_neighbor_affinity,
          R"pbdoc(Simulate a flockwork alpha-beta-model given an initial state as an edge
          list with varying reconnection rate alpha and varying disconnection rate beta
          (varying both over time and for each node). Choose neighbors of reacting nodes
          according to a weighted static network. Returns an instance of edge_changes.)pbdoc",
          py::arg("E"),
          py::arg("N"),
          py::arg("reconnection_rates"),
          py::arg("disconnection_rates"),
          py::arg("neighbor_affinity"),
          py::arg("t_run_total"),
          py::arg("tmax"),
          py::arg("use_random_rewiring") = false,
          py::arg("seed") = 0);

    m.def("equilibrate_flockwork_Q",
          &equilibrate_edgelist_seed,
          R"pbdoc(Equilibrates a flockwork Q-model given an initial state as an edge list.
          Returns new edge list.)pbdoc",
          py::arg("E"),
          py::arg("N"),
          py::arg("Q"),
          py::arg("seed"),
          py::arg("t_max") = 0,
          py::arg("use_Q_as_P") = false);

    m.def("equilibrate_flockwork_P",
          &equilibrate_edgelist_seed,
          R"pbdoc(Equilibrates a flockwork P-model given an initial state as an edge list.
          Returns new edge list.)pbdoc",
          py::arg("E"),
          py::arg("N"),
          py::arg("P"),
          py::arg("seed"),
          py::arg("t_max") = 0,
          py::arg("use_Q_as_P") = true);

    m.def("simulate_flockwork_Q",
          &simulate_flockwork,
          R"pbdoc(Simulates a flockwork Q-model given an initial state as an edge list.
          Returns new edge list.)pbdoc",
          py::arg("E"),
          py::arg("N"),
          py::arg("Q"),
          py::arg("seed"),
          py::arg("num_timesteps"));

    m.def("simulate_flockwork_P",
          &simulate_flockwork_P,
          R"pbdoc(Simulates a flockwork P-model given an initial state as an edge list.
          Returns new edge list.)pbdoc",
          py::arg("E"),
          py::arg("N"),
          py::arg("P"),
          py::arg("seed"),
          py::arg("num_timesteps"));

    m.def("simulate_flockwork_P_group_lifetimes",
          &simulate_flockwork_P_group_life_time,
          R"pbdoc(Simulates a flockwork P-model with initial state of all nodes
          unconnected. Returns a list of group life times.)pbdoc",
          py::arg("N"),
          py::arg("P"),
          py::arg("seed"),
          py::arg("num_timesteps"));

    m.def("measure_group_sizes_and_durations_for_edge_lists",
          &measure_group_sizes_and_durations,
          R"pbdoc(Provide an instance of tacoma.edge_lists and return a list of contact
          durations, a list of group size histograms and a list of durations lists,
          one for each group size.)pbdoc",
          py::arg("edge_lists"),
          py::arg("ignore_size_histograms") = false,
          py::arg("verbose") = false);

    m.def("measure_group_sizes_and_durations_for_edge_changes",
          &measure_group_sizes_and_durations_for_edge_changes,
          R"pbdoc(Provide an instance of tacoma.edge_lists and return a list of contact
          durations, a list of group size histograms and a list of durations
          lists, one for each group size.)pbdoc",
          py::arg("edge_changes"),
          py::arg("ignore_size_histogram_differences") = false,
          py::arg("verbose") = false);

    m.def("mean_degree_from_edge_lists",
          &mean_degree_from_edge_lists,
          R"pbdoc(Get a list of pairs where the first entry is a time point and the
          second entry is the mean degree of the current state, given a
          `_tacoma.edge_lists` instance.)pbdoc",
          py::arg("edge_lists"));

    m.def("mean_degree_from_edge_changes",
          &mean_degree_from_edge_changes,
          R"pbdoc(Get a list of pairs where the first entry
          is a time point and the second entry is the mean degree of the current
          state, given an `edge_changes` instance.)pbdoc",
          py::arg("edge_changes"));

    m.def("degree_distribution_from_edge_lists",
          &degree_distribution_from_edge_lists,
          R"pbdoc(Get a list of doubles where the k-th entry of the list is the
          time-averaged probability that a node has degree k.)pbdoc",
          py::arg("edge_lists"));

    m.def("degree_distribution_from_edge_changes",
          &degree_distribution_from_edge_changes,
          R"pbdoc(Get a list of doubles where the k-th entry of the list is the
          time-averaged probability that a node has degree k.)pbdoc",
          py::arg("edge_changes"));

    m.def("get_edge_counts",
          &get_edge_counts,
          R"pbdoc(Given an instance of :class:`_tacoma.edge_changes`, returns the lists
            `m_in`, `m_out` and `m`, which count the edge events as well as the
            edges at the given times. Note that `m` contains one entry more
            than `m_in` and `m_out` due to the initial edge list.)pbdoc",
          py::arg("edge_changes"));

    m.def("sample_from_edge_lists",
          &sample_from_edge_lists,
          R"pbdoc(Get a temporal network as a list of edge lists, given
            another instance of `edge_lists`, but resampled every dt.
            Alternatively, provide a number of time steps to divide (tmax-t0)
            into. if `sample_aggregates` is `True`, this does not sample the
            network state after each dt, but rather gives a graph of all edges
            being present in the last time step. However, you should use
            :func:`_tacoma.bin_from_edge_lists` instead.)pbdoc",
          py::arg("edge_lists"),
          py::arg("dt") = 0.0,
          py::arg("N_time_steps") = 0,
          py::arg("sample_aggregates") = false,
          py::arg("verbose") = false);

    m.def("sample_from_edge_changes",
          &sample_from_edge_changes,
          R"pbdoc(Get a temporal network as a :class:`_tacoma.edge_lists`,
            given an instance of :class:`_tacoma.edge_changes`,
            but resampled every dt. Alternatively, provide a number of time steps
            to divide (tmax-t0) into. if `sample_aggregates` is `True`, this does
            not sample the network state after each dt, but rather gives a graph
            of all edges being present in the last time step. However, you should
            use :func:`_tacoma.bin_from_edge_changes instead`.)pbdoc",
          py::arg("edge_changes"),
          py::arg("dt") = 0.0,
          py::arg("N_time_steps") = 0,
          py::arg("sample_aggregates") = false,
          py::arg("verbose") = false);

    m.def("bin_from_edge_lists",
          &bin_from_edge_lists,
          R"pbdoc(Get a temporal network as a list of edge lists, given
            another instance of `edge_lists`, but binned for every dt. Alternatively,
            provide a number of time steps to divide (tmax-t0) into. This
            does not sample the network state after each dt, but rather gives a
            graph of all edges being present in the last time step.)pbdoc",
          py::arg("edge_lists"),
          py::arg("dt") = 0.0,
          py::arg("N_time_steps") = 0,
          py::arg("verbose") = false);

    m.def("bin_from_edge_changes",
          &bin_from_edge_changes,
          R"pbdoc(Get a temporal network as :class:`_tacoma.edge_lists`,
            given an instance of `edge_changes`, but binned every dt.
            Alternatively, provide a number of time steps to divide (tmax-t0) into.
            This does not sample the network state after each dt, but rather gives a
            graph of all edges being present in the last time step.)pbdoc",
          py::arg("edge_changes"),
          py::arg("dt") = 0.0,
          py::arg("N_time_steps") = 0,
          py::arg("verbose") = false);

    m.def("slice_edge_lists",
          &slice_edge_lists,
          R"pbdoc(Get a temporal network as a list of edge lists,
            given another instance of :class:`_tacoma.edge_lists`,
            but only for times ``new_t0 <= t < new_tmax``.)pbdoc",
          py::arg("edge_lists"),
          py::arg("new_t0"),
          py::arg("new_tmax"),
          py::arg("verbose") = false);

    m.def("slice_edge_changes",
          &slice_edge_changes,
          R"pbdoc(Get a temporal network as a list of edge changes,
            given another instance of :class:`_tacoma.edge_changes`,
            but only for times ``new_t0 <= t < new_tmax``.)pbdoc",
          py::arg("edge_changes"),
          py::arg("new_t0"),
          py::arg("new_tmax"),
          py::arg("verbose") = false);

    m.def("binned_social_trajectory_from_edge_changes",
          &binned_social_trajectory_from_edge_changes,
          R"pbdoc(Get a social trajectory of a node, binned every dt
            (list of group indices this node was part of
            in :math:`[t,t+\Delta t)` ).)pbdoc",
          py::arg("edge_changes"),
          py::arg("node"),
          py::arg("dt") = 0.0,
          py::arg("N_time_steps") = 0,
          py::arg("verbose") = false);

    m.def("social_trajectory_from_edge_changes",
          &social_trajectory_from_edge_changes,
          "Get a social trajectory of a node.",
          py::arg("edge_changes"),
          py::arg("node"),
          py::arg("verbose") = false);

    m.def("social_trajectory_from_edge_lists",
          &social_trajectory_from_edge_lists,
          "Get a social trajectory of a node.",
          py::arg("edge_lists"),
          py::arg("node"),
          py::arg("verbose") = false);

    m.def("edge_trajectories_from_edge_lists",
          &edge_trajectories_from_edge_lists,
          "For each edge, get all time pairs where the edge existed.",
          py::arg("edge_lists"),
          py::arg("return_edge_similarities") = false,
          py::arg("verbose") = false);

    m.def("edge_trajectories_from_edge_changes",
          &edge_trajectories_from_edge_changes,
          R"pbdoc(For each edge, get all time pairs where the
            edge existed, given an instance of `edge_changes`.)pbdoc",
          py::arg("edge_changes"),
          py::arg("return_edge_similarities") = false,
          py::arg("verbose") = false);

    m.def("convert_edge_trajectories_to_edge_changes",
          &edge_trajectories_to_edge_changes,
          R"pbdoc(Convert a :class:`_tacoma.edge_trajectories`
            instance to an instance of :class:`_tacoma.edge_changes`.)pbdoc",
          py::arg("traj"));

    m.def("verify_edge_lists",
          &verify_edge_lists,
          R"pbdoc(Verify that the given temporal network is in a valid format.
            Returns the number of rule violations and is verbose about
            what's wrong.)pbdoc",
          py::arg("edge_lists"),
          py::arg("verbose") = false);

    m.def("verify_edge_changes",
          &verify_edge_changes,
          R"pbdoc(Verify that the given temporal network is in a valid format.
            Returns the number of rule violations and is verbose about
            what's wrong.)pbdoc",
          py::arg("edge_changes"),
          py::arg("verbose") = false);

    m.def("convert_edge_lists",
          &convert_edge_lists,
          R"pbdoc(Convert an instance of :class:`_tacoma.edge_lists`
            to an instance of :class:`_tacoma.edge_changes`.)pbdoc",
          py::arg("edge_lists"),
          py::arg("verbose") = false);

    m.def("convert_edge_changes",
          &convert_edge_changes,
          R"pbdoc(Convert an instance of :class:`_tacoma.edge_changes` to an
            instance of :class:`_tacoma.edge_lists`.)pbdoc",
          py::arg("edge_changes"),
          py::arg("verbose") = false);

    m.def("concatenate_edge_lists",
          &concatenate_edge_lists,
          R"pbdoc(Concatenate a list of :class:`_tacoma.edge_lists` to a
            single instance of :class:`_tacoma.edge_lists`.)pbdoc",
          py::arg("list_of_edge_lists"),
          py::arg("verbose") = false);

    m.def("concatenate_edge_changes",
          &concatenate_edge_changes,
          R"pbdoc(Concatenate a list of :class:`_tacoma.edge_changes`
            to a single instance of :class:`_tacoma.edge_changes`.)pbdoc",
          py::arg("list_of_edge_changes"),
          py::arg("verbose") = false);

    m.def("binned_social_trajectory_from_edge_lists",
          &binned_social_trajectory_from_edge_lists,
          R"pbdoc(Get a social trajectory of a node, given an instance of
            :class:`_tacoma.edge_lists`, binned every dt (list of group indices
            this node was part of in :math:`[t,t+dt)`).)pbdoc",
          py::arg("edge_lists"),
          py::arg("node"),
          py::arg("dt") = 0.0,
          py::arg("N_time_steps") = 0,
          py::arg("verbose") = false);

    m.def("get_flockwork_P_args",
          &get_flockwork_P_args,
          R"pbdoc(Calculate the rewiring_rate :math:`\gamma(t)` and probability to
            reconnect :math:`P(t)` as well as the other important parameters to
            simulate a flockwork P-model with varying rates.)pbdoc",
          py::arg("edge_changes"),
          py::arg("dt") = 0.0,
          py::arg("N_time_steps") = 0,
          py::arg("k_over_k_real_scaling") = 1.0,
          py::arg("gamma_scaling") = 1.0,
          py::arg("P_scaling") = 1.0,
          py::arg("aggregated_network") = map<pair<size_t, size_t>, double>(),
          py::arg("ensure_empty_network") = false,
          py::arg("adjust_last_bin_if_dt_does_not_fit") = false,
          py::arg("verbose") = false);

    m.def("get_flockwork_alpha_beta_args", &get_flockwork_alpha_beta_args,
          R"pbdoc(Calculate the reconnection rate :math:`\alpha(t)` and disconnection rate
            :math:`\beta(t)` as well as the other important parameters to simulate a
            flockwork-alpha-beta-model with varying rates.)pbdoc",
          py::arg("edge_changes"),
          py::arg("dt") = 0.0,
          py::arg("N_time_steps") = 0,
          py::arg("k_over_k_real_scaling") = 1.0,
          py::arg("alpha_scaling") = 1.0,
          py::arg("beta_scaling") = 1.0,
          py::arg("aggregated_network") = map<pair<size_t, size_t>, double>(),
          py::arg("ensure_empty_network") = false,
          py::arg("consider_looped_network_to_get_final_events") = false,
          py::arg("verbose") = false);

    m.def("get_flockwork_P_node_parameters_gamma_and_P", get_node_gamma_and_P,
          R"pbdoc(Calculate the mean node specific activity :math:`\gamma_i`
            and reconnection probability :math:`P_i` for each node.)pbdoc",
          py::arg("edge_changes"),
          py::arg("gamma"),
          py::arg("P"),
          py::arg("use_event_rate_method") = false);

    m.def("get_flockwork_node_parameters_alpha_and_beta", get_node_alpha_and_beta,
          R"pbdoc(Calculate the mean node specific reconnection rate factor
            :math:`\alpha_i` and disconnection rate factor :math`\beta_i` for each node.)pbdoc",
          py::arg("edge_changes"),
          py::arg("reconnection_rate"),
          py::arg("disconnecton_rate"),
          py::arg("k_over_k_real_scaling") = 1.0,
          py::arg("apply_mean_correction") = true,
          py::arg("verbose") = false);

    m.def("get_flockwork_node_rates_alpha_and_beta", get_time_dependent_node_alpha_and_beta,
          R"pbdoc(Calculate the mean node specific reconnection rate factor
            :math:`\alpha_i` and disconnection rate factor :math:`\beta_i` for each node.)pbdoc",
          py::arg("edge_changes"),
          py::arg("reconnection_rate"),
          py::arg("disconnecton_rate"),
          py::arg("apply_mean_correction") = true);

    m.def("get_flockwork_P_node_parameters_alpha_and_beta_from_gamma_and_P", get_node_gamma_and_P,
          R"pbdoc(Calculate the mean node specific activity :math:`\gamma_i` and reconnection
            probability :math:`P_i` for each node.)pbdoc",
          py::arg("edge_changes"),
          py::arg("gamma"),
          py::arg("P"),
          py::arg("use_event_rate_method") = true);

    m.def("estimate_k_scaling_gradient_descent", &estimate_k_scaling_gradient_descent,
          R"pbdoc(Estimate the scaling of $\left\langle k\right\rangle$ that's
            necessary to revert the effects of binning using gradient descent.)pbdoc",
          py::arg("edge_changes"),
          py::arg("dt_for_inference"),
          py::arg("dt_for_binning"),
          py::arg("measurements_per_configuration") = 6,
          py::arg("learning_rate") = 0.5,
          py::arg("relative_error") = 1e-2,
          py::arg("N_eval_max") = 100,
          py::arg("verbose") = true);

    m.def("estimate_k_scaling_gradient_descent_RMSE", &estimate_k_scaling_gradient_descent_RMSE,
          R"pbdoc(Estimate the scaling of $\left\langle k\right\rangle$ that's
            necessary to revert the effects of binning using gradient descent.)pbdoc",
          py::arg("edge_changes"),
          py::arg("dt_for_inference"),
          py::arg("dt_for_binning"),
          py::arg("measurements_per_configuration") = 6,
          py::arg("learning_rate") = 0.5,
          py::arg("relative_error") = 1e-2,
          py::arg("N_eval_max") = 100,
          py::arg("verbose") = true);

    m.def("ZSBB_model", &ZSBB_model,
          "Simulate model after Zhao, Stehle, Bianconi, and Barrat.",
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
          py::arg("verbose") = false);

    m.def("dynamic_RGG", &dynamic_RGG,
          "Simulate a dynamic random geometric graph model.",
          py::arg("N"),
          py::arg("t_run_total"),
          py::arg("step_distance") = 0.0,
          py::arg("mean_link_duration") = 0.0,
          py::arg("critical_density") = 0.65,
          py::arg("periodic_boundary_conditions_for_link_building") = false,
          py::arg("record_sizes_and_durations") = false,
          py::arg("seed") = 0,
          py::arg("verbose") = false);

    m.def("gillespie_SIS_on_edge_lists", &gillespie_on_edge_lists<SIS>,
          "Perform a Gillespie SIS simulation on edge lists.",
          py::arg("edge_lists"),
          py::arg("SIS"),
          py::arg("is_static") = false,
          py::arg("verbose") = false);

    m.def("gillespie_SIS_on_edge_changes", &gillespie_on_edge_changes<SIS>,
          "Perform a Gillespie SIS simulation on edge changes.",
          py::arg("edge_changes"),
          py::arg("SIS"),
          py::arg("verbose") = false);

    m.def("markov_SIS_on_edge_lists", &markov_on_edge_lists<MARKOV_SIS>,
          "Perform a Markov SIS simulation on edge lists.",
          py::arg("edge_lists"),
          py::arg("MARKOV_SIS"),
          py::arg("dt"),
          py::arg("is_static") = false,
          py::arg("verbose") = false);

    m.def("markov_SIS_on_edge_changes", &markov_on_edge_changes<MARKOV_SIS>,
          "Perform a Markov SIS simulation on edge changes.",
          py::arg("edge_changes"),
          py::arg("MARKOV_SIS"),
          py::arg("dt"),
          py::arg("verbose") = false);

    m.def("gillespie_QS_SIS_on_edge_lists", &gillespie_on_edge_lists<QS_SIS>,
          "Perform a quasi-stationary Gillespie SIS simulation on edge lists.",
          py::arg("edge_lists"),
          py::arg("QS_SIS"),
          py::arg("is_static") = false,
          py::arg("verbose") = false);

    m.def("gillespie_QS_SIS_on_edge_changes", &gillespie_on_edge_changes<QS_SIS>,
          "Perform a quasi-stationary Gillespie SIS simulation on edge changes.",
          py::arg("edge_changes"),
          py::arg("QS_SIS"),
          py::arg("verbose") = false);

    m.def("gillespie_eSIS_on_edge_lists", &gillespie_on_edge_lists<eSIS>,
          R"pbdoc(Perform a Gillespie :math:`\varepsilon`-SIS simulation on edge lists.)pbdoc",
          py::arg("edge_lists"),
          py::arg("eSIS"),
          py::arg("is_static") = false,
          py::arg("verbose") = false);

    m.def("gillespie_eSIS_on_edge_changes", &gillespie_on_edge_changes<eSIS>,
          R"pbdoc(Perform a Gillespie :math:`\varepsilon`-SIS simulation on edge changes.)pbdoc",
          py::arg("edge_changes"),
          py::arg("eSIS"),
          py::arg("verbose") = false);

    m.def("gillespie_coverage_SIS_on_edge_changes", &gillespie_on_edge_changes<coverage_SIS>,
          R"pbdoc(Perform a Gillespie coverage-SIS simulation on edge changes.)pbdoc",
          py::arg("edge_changes"),
          py::arg("coverage_SIS"),
          py::arg("verbose") = false);

    m.def("gillespie_coverage_SIS_on_edge_lists", &gillespie_on_edge_lists<coverage_SIS>,
          R"pbdoc(Perform a Gillespie coverage-SIS simulation on edge lists.)pbdoc",
          py::arg("edge_lists"),
          py::arg("coverage_SIS"),
          py::arg("is_static") = false,
          py::arg("verbose") = false);

    m.def("gillespie_cluster_size_SIS_on_edge_changes", &gillespie_on_edge_changes<cluster_size_SIS>,
          R"pbdoc(Perform a Gillespie cluster-size-SIS simulation on edge changes.)pbdoc",
          py::arg("edge_changes"),
          py::arg("cluster_size_SIS"),
          py::arg("verbose") = false);

    m.def("gillespie_cluster_size_SIS_on_edge_lists", &gillespie_on_edge_lists<cluster_size_SIS>,
          R"pbdoc(Perform a Gillespie cluster-size-SIS simulation on edge lists.)pbdoc",
          py::arg("edge_lists"),
          py::arg("cluster_size_SIS"),
          py::arg("is_static") = false,
          py::arg("verbose") = false);


    m.def("gillespie_node_based_SIS_on_edge_changes", &gillespie_on_edge_changes<node_based_SIS>,
          "Perform a Gillespie SIS simulation on edge changes.",
          py::arg("edge_changes"),
          py::arg("SIS"),
          py::arg("verbose") = false);

    m.def("gillespie_node_based_SIS_on_edge_lists", &gillespie_on_edge_lists<node_based_SIS>,
          "Perform a Gillespie SIS simulation on edge lists.",
          py::arg("edge_lists"),
          py::arg("SIS"),
          py::arg("is_static") = false,
          py::arg("verbose") = false);

    m.def("gillespie_SI_on_edge_lists", &gillespie_on_edge_lists<SI>,
          "Perform a Gillespie SI simulation on edge lists.",
          py::arg("edge_lists"),
          py::arg("SI"),
          py::arg("is_static") = false,
          py::arg("verbose") = false);

    m.def("gillespie_SI_on_edge_changes", &gillespie_on_edge_changes<SI>,
          "Perform a Gillespie SI simulation on edge changes.",
          py::arg("edge_changes"),
          py::arg("SI"),
          py::arg("verbose") = false);

    m.def("gillespie_SIR_on_edge_lists", &gillespie_on_edge_lists<SIR>,
          "Perform a Gillespie SIR simulation on edge lists.",
          py::arg("edge_lists"),
          py::arg("SIR"),
          py::arg("is_static") = false,
          py::arg("verbose") = false);

    m.def("gillespie_SIR_on_edge_changes", &gillespie_on_edge_changes<SIR>,
          "Perform a Gillespie SIR simulation on edge changes.",
          py::arg("edge_changes"),
          py::arg("SIR"),
          py::arg("verbose") = false);

    m.def("gillespie_SIRS_on_edge_lists", &gillespie_on_edge_lists<SIRS>,
          "Perform a Gillespie SIRS simulation on edge lists.",
          py::arg("edge_lists"),
          py::arg("SIRS"),
          py::arg("is_static") = false,
          py::arg("verbose") = false);

    m.def("gillespie_SIRS_on_edge_changes", &gillespie_on_edge_changes<SIRS>,
          "Perform a Gillespie SIRS simulation on edge changes.",
          py::arg("edge_changes"),
          py::arg("SIRS"),
          py::arg("verbose") = false);

    m.def("gillespie_SIS_on_EdgeActivityModel", &gillespie_on_model<EdgeActivityModel,SIS>,
          "Perform a Gillespie SIS simulation on the edge activity model.",
          py::arg("edge_activity_model"),
          py::arg("SIS"),
          py::arg("reset_simulation_objects") = true,
          py::arg("verbose") = false);

    m.def("simulate_EdgeActivityModel", &gillespie_model_simulation<EdgeActivityModel>,
          "Perform a Gillespie simulation of the edge activity model.",
          py::arg("model"),
          py::arg("t_simulation"),
          py::arg("reset_simulation_objects") = true,
          py::arg("verbose") = false);

    m.def("simulate_FlockworkPModel", &gillespie_model_simulation<FlockworkPModel>,
          "Perform a Gillespie simulation of the Flockwork P-model.",
          py::arg("model"),
          py::arg("t_simulation"),
          py::arg("reset_simulation_objects") = true,
          py::arg("verbose") = false);

    m.def("gillespie_coverage_SIS_on_EdgeActivityModel", &gillespie_on_model<EdgeActivityModel,coverage_SIS>,
          "Perform a Gillespie coverage SIS simulation on the edge activity model.",
          py::arg("edge_activity_model"),
          py::arg("coverage_SIS"),
          py::arg("reset_simulation_objects") = true,
          py::arg("verbose") = false);

    m.def("gillespie_cluster_size_SIS_on_EdgeActivityModel", &gillespie_on_model<EdgeActivityModel,cluster_size_SIS>,
          "Perform a Gillespie cluster size SIS simulation on the edge activity model.",
          py::arg("edge_activity_model"),
          py::arg("cluster_size_SIS"),
          py::arg("reset_simulation_objects") = true,
          py::arg("verbose") = false);

    m.def("gillespie_QS_SIS_on_EdgeActivityModel", &gillespie_on_model<EdgeActivityModel,QS_SIS>,
          "Perform a quasi-stationary Gillespie SIS simulation on the edge activity model.",
          py::arg("edge_activity_model"),
          py::arg("QS_SIS"),
          py::arg("reset_simulation_objects") = false,
          py::arg("verbose") = false);

    m.def("gillespie_eSIS_on_EdgeActivityModel", &gillespie_on_model<EdgeActivityModel,eSIS>,
          R"pbdoc(Perform a Gillespie :math:`\varepsilon`-SIS simulation on the edge activity model.)pbdoc",
          py::arg("activity_model"),
          py::arg("eSIS"),
          py::arg("reset_simulation_objects") = true,
          py::arg("verbose") = false);

    m.def("gillespie_SIR_on_EdgeActivityModel", &gillespie_on_model<EdgeActivityModel,SIR>,
          "Perform a Gillespie SIR simulation on the edge activity model.",
          py::arg("edge_activity_model"),
          py::arg("SIR"),
          py::arg("reset_simulation_objects") = true,
          py::arg("verbose") = false);

    m.def("gillespie_SIRS_on_EdgeActivityModel", &gillespie_on_model<EdgeActivityModel,SIRS>,
          "Perform a Gillespie SIRS simulation on the edge activity model.",
          py::arg("edge_activity_model"),
          py::arg("SIRS"),
          py::arg("reset_simulation_objects") = true,
          py::arg("verbose") = false);

    m.def("gillespie_SI_on_EdgeActivityModel", &gillespie_on_model<EdgeActivityModel,SI>,
          "Perform a Gillespie SI simulation on the edge activity model.",
          py::arg("edge_activity_model"),
          py::arg("SI"),
          py::arg("reset_simulation_objects") = true,
          py::arg("verbose") = false);

    m.def("gillespie_node_based_SIS_on_EdgeActivityModel", &gillespie_on_model<EdgeActivityModel,node_based_SIS>,
          "Perform a node-based Gillespie SIS simulation on the edge activity model.",
          py::arg("edge_activity_model"),
          py::arg("SI"),
          py::arg("reset_simulation_objects") = true,
          py::arg("verbose") = false);

    m.def("markov_SIS_on_EdgeActivityModel", &markov_on_model<EdgeActivityModel,MARKOV_SIS>,
          "Perform a mixed Markov-Gillespie SIS simulation on the edge activity model.",
          py::arg("edge_activity_model"),
          py::arg("MARKOV_SIS"),
          py::arg("max_dt"),
          py::arg("reset_simulation_objects") = true,
          py::arg("verbose") = false);

    m.def("gillespie_SIS_on_FlockworkPModel", &gillespie_on_model<FlockworkPModel,SIS>,
          "Perform a Gillespie SIS simulation on the Flockwork P-model.",
          py::arg("model"),
          py::arg("SIS"),
          py::arg("reset_simulation_objects") = true,
          py::arg("verbose") = false);

    m.def("gillespie_coverage_SIS_on_FlockworkPModel", &gillespie_on_model<FlockworkPModel,coverage_SIS>,
          "Perform a Gillespie coverage-SIS simulation on the Flockwork P-model.",
          py::arg("model"),
          py::arg("coverage_SIS"),
          py::arg("reset_simulation_objects") = true,
          py::arg("verbose") = false);

    m.def("gillespie_cluster_size_SIS_on_FlockworkPModel", &gillespie_on_model<FlockworkPModel,cluster_size_SIS>,
          "Perform a Gillespie cluster-size-SIS simulation on the Flockwork P-model.",
          py::arg("model"),
          py::arg("cluster_size_SIS"),
          py::arg("reset_simulation_objects") = true,
          py::arg("verbose") = false);

    m.def("gillespie_QS_SIS_on_FlockworkPModel", &gillespie_on_model<FlockworkPModel,QS_SIS>,
          "Perform a quasi-stationary Gillespie SIS simulation on the Flockwork P-model.",
          py::arg("model"),
          py::arg("QS_SIS"),
          py::arg("reset_simulation_objects") = false,
          py::arg("verbose") = false);

    m.def("gillespie_eSIS_on_FlockworkPModel", &gillespie_on_model<FlockworkPModel,eSIS>,
          R"pbdoc(Perform a Gillespie :math:`\varepsilon`-SIS simulation on the Flockwork P-model.)pbdoc",
          py::arg("model"),
          py::arg("eSIS"),
          py::arg("reset_simulation_objects") = true,
          py::arg("verbose") = false);

    m.def("gillespie_SIR_on_FlockworkPModel", &gillespie_on_model<FlockworkPModel,SIR>,
          "Perform a Gillespie SIR simulation on the Flockwork P-model.",
          py::arg("model"),
          py::arg("SIR"),
          py::arg("reset_simulation_objects") = true,
          py::arg("verbose") = false);

    m.def("gillespie_SIRS_on_FlockworkPModel", &gillespie_on_model<FlockworkPModel,SIRS>,
          "Perform a Gillespie SIRS simulation on the Flockwork P-model.",
          py::arg("model"),
          py::arg("SIRS"),
          py::arg("reset_simulation_objects") = true,
          py::arg("verbose") = false);

    m.def("gillespie_SI_on_FlockworkPModel", &gillespie_on_model<FlockworkPModel,SI>,
          "Perform a Gillespie SI simulation on the Flockwork P-model.",
          py::arg("model"),
          py::arg("SI"),
          py::arg("reset_simulation_objects") = true,
          py::arg("verbose") = false);

    m.def("gillespie_node_based_SIS_on_FlockworkPModel", &gillespie_on_model<FlockworkPModel,node_based_SIS>,
          "Perform a node-based Gillespie SIS simulation on the Flockwork P-model.",
          py::arg("model"),
          py::arg("SI"),
          py::arg("reset_simulation_objects") = true,
          py::arg("verbose") = false);

    m.def("markov_SIS_on_FlockworkPModel", &markov_on_model<FlockworkPModel,MARKOV_SIS>,
          "Perform a mixed Markov-Gillespie SIS simulation on the Flockwork P-model.",
          py::arg("model"),
          py::arg("MARKOV_SIS"),
          py::arg("max_dt"),
          py::arg("reset_simulation_objects") = true,
          py::arg("verbose") = false);    


    py::class_<edge_changes>(m, "edge_changes", R"pbdoc(Description of a temporal network by listing the changes of edges at a certain time.)pbdoc")
        .def(py::init<>())
        .def(py::init<const edge_changes_with_histograms &>(),
             py::arg("edge_changes_with_histograms"),
             "Initialize from an instance of :class:`_tacoma.edge_changes_with_histograms`")
        .def("copy_from", &edge_changes::copy_from, R"pbdoc([deprecated] copy the relevant information from an instance of :class:`_tacoma.edge_changes_with_histograms`)pbdoc")
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

    py::class_<edge_lists>(m, "edge_lists", R"pbdoc(
            Description of a temporal network by listing
            the edges of the networks at a certain time.
            )pbdoc")
        .def(py::init<>())
        .def(py::init<const edge_lists_with_histograms &>(),
             py::arg("edge_lists_with_histograms"),
             "Initialize from an instance of :class:`_tacoma.edge_lists_with_histograms`")
        .def("copy_from", &edge_lists::copy_from, R"pbdoc([deprecated] copy the relevant information from an instance of :class:`_tacoma.edge_lists_with_histograms`)pbdoc")
        .def_readwrite("int_to_node", &edge_lists::int_to_node, R"pbdoc(A dictionary int -> string which keeps the original node names.)pbdoc")
        .def_readwrite("time_unit", &edge_lists::time_unit, R"pbdoc(A string containing the unit of time for this network.)pbdoc")
        .def_readwrite("notes", &edge_lists::notes, R"pbdoc(A string containing additional notes for this network.)pbdoc")
        .def_readwrite("t", &edge_lists::t, R"pbdoc(An ordered list containing the time points at which the new edge list becomes active.)pbdoc")
        .def_readwrite("edges", &edge_lists::edges, R"pbdoc(A list containing the edge list of the network at the correspoding time in `t`)pbdoc")
        .def_readwrite("N", &edge_lists::N, R"pbdoc(Number of nodes)pbdoc")
        .def_readwrite("tmax", &edge_lists::tmax, R"pbdoc(The final time)pbdoc");

    py::class_<edge_lists_with_histograms>(m, "edge_lists_with_histograms",
                                           R"pbdoc(
                Similar to the :class:`_tacoma.edge_lists` class but with additional analysis results.
            )pbdoc")
        .def(py::init<>())
        .def_readwrite("t", &edge_lists_with_histograms::t,
                       R"pbdoc(
                An ordered list containing the time points at which the new edge list becomes active.
            )pbdoc")
        .def_readwrite("edges", &edge_lists_with_histograms::edges,
                       R"pbdoc(
                A list containing the edge list of the network at the correspoding time in `t`
            )pbdoc")
        .def_readwrite("size_histograms", &edge_lists_with_histograms::size_histograms,
                       R"pbdoc(
                A list of :obj:`dict` s, one dictionary for each edge list, containing the number of observed groups of size `m` at this time.
            )pbdoc")
        .def_readwrite("group_durations", &edge_lists_with_histograms::group_durations,
                       R"pbdoc(
                A list of :obj:`list` s where the `m`-th entry contains all durations of all observed groups of size `m`.
            )pbdoc")
        .def_readwrite("N", &edge_lists_with_histograms::N, "The number of nodes.")
        .def_readwrite("tmax", &edge_lists_with_histograms::tmax, "The maximum time.");

    py::class_<edge_changes_with_histograms>(m, "edge_changes_with_histograms",
                                             R"pbdoc(
                Similar to the :class:`_tacoma.edge_changes` class but with additional analysis results.
            )pbdoc")
        .def(py::init<>())
        .def_readwrite("t", &edge_changes_with_histograms::t, R"pbdoc(An ordered list containing the time points at which changes occur.)pbdoc")
        .def_readwrite("edges_in", &edge_changes_with_histograms::edges_in, R"pbdoc(A list containing the edges coming in to the network at the correspoding time in `t`.)pbdoc")
        .def_readwrite("edges_out", &edge_changes_with_histograms::edges_out, R"pbdoc(A list containing the edges leaving the network at the correspoding time in `t`.)pbdoc")
        .def_readwrite("initial_size_histogram", &edge_changes_with_histograms::initial_size_histogram,
                       R"pbdoc(
                A :obj:`dict` containing group sizes as keys and counts of groups of corresponding sizes as values.
            )pbdoc")
        .def_readwrite("group_changes", &edge_changes_with_histograms::group_changes,
                       R"pbdoc(
                A list :obj:`dict` , one dict for each time step, containing group sizes as keys and change of value to the previous time step as values.
            )pbdoc")
        .def_readwrite("final_size_histogram", &edge_changes_with_histograms::final_size_histogram,
                       R"pbdoc(
                A :obj:`dict` containing group sizes as keys and counts of groups of corresponding sizes as values.
            )pbdoc")
        .def_readwrite("contact_durations", &edge_changes_with_histograms::contact_durations,
                       R"pbdoc(
                A list of :obj:`int` s, each entry is the duration of a single contact (in time steps).
            )pbdoc")
        .def_readwrite("inter_contact_durations", &edge_changes_with_histograms::inter_contact_durations,
                       R"pbdoc(
                A list of :obj:`int` s, each entry is the number of time steps it took a single node to reconnect.
            )pb
            doc")
        .def_readwrite("group_durations", &edge_changes_with_histograms::group_durations,
            R"pbdoc(
                A list of :obj:`list` s, the `m`-th entry contains all durations of all observed groups of size `m`.
            )pbdoc")
        .def_readwrite("N", &edge_changes_with_histograms::N, "The number of nodes.")
        .def_readwrite("t0", &edge_changes_with_histograms::t0, "The initial time.")
        .def_readwrite("tmax", &edge_changes_with_histograms::tmax, "The final time.")
        .def_readwrite("edges_initial", &edge_changes_with_histograms::edges_initial, R"pbdoc(A list containing the edges of the network at time `t0`)pbdoc");

    py::class_<group_sizes_and_durations>(m, "group_sizes_and_durations")
        .def(py::init<>())
        .def_readwrite("contact_durations", &group_sizes_and_durations::contact_durations,
                       R"pbdoc(
                A list of :obj:`float` s, each entry is the duration of a single contact.
            )pbdoc")
        .def_readwrite("size_histograms", &group_sizes_and_durations::size_histograms,
                       R"pbdoc(
                A list of :obj:`dict` s, one dictionary for each edge list, containing the number of observed groups of size `m` at this time. Has only one entry if an instance of tacoma.edge_changes was provided.
            )pbdoc")
        .def_readwrite("size_histogram_differences", &group_sizes_and_durations::size_histogram_differences,
                       R"pbdoc(
                A list :obj:`dict` , one dict for each time step, containing
                group sizes as keys and change of value to the previous time step as values.
            )pbdoc")
        .def_readwrite("group_durations", &group_sizes_and_durations::group_durations,
                       R"pbdoc(
                A list of :obj:`list` s, the :math:`m`-th entry contains
                all durations of all observed groups of size $m$.
            )pbdoc")
        .def_readwrite("aggregated_size_histogram", &group_sizes_and_durations::aggregated_size_histogram,
                       R"pbdoc(
                A list of :obj:`list` s, the :math:`m`-th entry contains
                the time-averaged number of groups of size :math:`m`,
                :math:`\overline{N_m}=(1/t_\mathrm{max})\int_0^{t_\mathrm{max}}dt N_m(t)`.
            )pbdoc")
        .def_readwrite("aggregated_network", &group_sizes_and_durations::aggregated_network,
                       R"pbdoc(
                A :obj:`dict` where each (key, value)-pair is an edge
                and its corresponding total time it was active.
            )pbdoc");

    py::class_<edge_weight>(m, "edge_weight",
                            R"pbdoc(Helper class for internal usage. Creates `value = 0`
            when initiated such that one can easily use it as a counter
            in a map without checking whether or not the object exists.)pbdoc")
        .def(py::init<>())
        .def_readwrite("value", &edge_weight::value);

    /*
    py::class_<EPI>(m, "EPI",
                            R"pbdoc(Helper class for internal usage. Maps the epidemic state of
                            a node to an integer.)pbdoc")
        .def(py::init<>())
        .def_read("S", &EPI::S, "A node is susceptible to disease.")
        .def_read("I", &EPI::I, "A node is infectious.")
        .def_read("R", &EPI::R, "A node is recovered.")
        .def_read("E", &EPI::E, "A node is exposed (infected but not infectious).")
        .def_read("V", &EPI::V, R"pbdoc(A node is vaccinated (This means the same as ``R`` in 
                                        most models, but in SIRS it implies that vaccinated nodes 
                                        cannot have waning immunity).)pbdoc");
    */

    py::class_<social_trajectory_entry>(m, "social_trajectory_entry",
                                        "Each :class:`_tacoma.social_trajectory_entry` ")
        .def(py::init<>())
        .def_readwrite("hash", &social_trajectory_entry::hash,
                       "A group identifier which has low probability of doubling.")
        .def_readwrite("size", &social_trajectory_entry::size,
                       "Number of nodes within this group.")
        .def_readwrite("time_pairs", &social_trajectory_entry::time_pairs,
                       R"pbdoc(List[Tuple[double, double]] ordered time pairs denoting
                       the time intervals in which the group existed)pbdoc");

    py::class_<edge_trajectory_entry>(m, "edge_trajectory_entry",
                                      R"pbdoc(This is an entry of an edge-based notation of a temporal
        network. Instead of getting lists of edges
        ordered in time, a list of :class:`_tacoma.edge_trajectory_entry` contains,
        for each edge, a list of time pairs denoting the times the edge exists.)pbdoc")
        .def(py::init<>())
        .def(py::init<
                 const pair<size_t, size_t> &,
                 const vector<pair<double, double>> &>(),
             py::arg("edge"),
             py::arg("time_pairs"),
             "Initialize from a pair of ints and a list of pairs of doubles.")
        .def_readwrite("edge", &edge_trajectory_entry::edge, R"pbdoc(Tuple[int, int] containing the nodes belonging to this edge.)pbdoc")
        .def_readwrite("time_pairs", &edge_trajectory_entry::time_pairs, R"pbdoc(
            List[Tuple[double, double]]
            A list containing ordered time pairs :math:`(t_i^{(i)}, t_i^{(f)})`, where each time pair contains
            the time point :math:`t_i^{(i)}` when the edge is switched on (created) and
            the time :math:`t_i^{(f)}` when the edge is switched off (deleted).
            )pbdoc");

    py::class_<flockwork_args>(m, "flockwork_args", py::dynamic_attr(), R"pbdoc(
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
        .def_readwrite("new_time", &flockwork_args::new_time, "The bin edges of the new time bins");

    py::class_<flockwork_alpha_beta_args>(m, "flockwork_alpha_beta_args", py::dynamic_attr(), R"pbdoc(
                    The arguments which are passed to the flockwork alpha-beta simulation
                    function. An instance of this is returned by the parameter estimation procedure.
                )pbdoc")
        .def(py::init<>())
        .def_readwrite("E", &flockwork_alpha_beta_args::E, "Edge list of the initial state")
        .def_readwrite("N", &flockwork_alpha_beta_args::N, "Number of nodes")
        .def_readwrite("disconnection_rate", &flockwork_alpha_beta_args::disconnection_rate, R"pbdoc(A list of floats describing the time-dependent disconnection rate (corresponding times in `reconnection_rate`))pbdoc")
        .def_readwrite("reconnection_rate", &flockwork_alpha_beta_args::reconnection_rate, R"pbdoc(A list of pairs of doubles, each entry contains the time and the reconnection rate per node)pbdoc")
        .def_readwrite("neighbor_affinity", &flockwork_alpha_beta_args::neighbor_affinity, R"pbdoc(
        A list, for each node contains two lists, one containing the node's
        neighbors and the second one containing the node affinity between them.)pbdoc")
        .def_readwrite("tmax", &flockwork_alpha_beta_args::tmax, R"pbdoc(The time at which the last value of `P` and the `rewiring_rate` changes (i.e. the maximum time until the parameteres are defined))pbdoc")
        .def_readwrite("M", &flockwork_alpha_beta_args::M, "The integral over number of edges in the last time interval.")
        .def_readwrite("m", &flockwork_alpha_beta_args::m, "The number of edges in the last time interval.")
        .def_readwrite("m_in", &flockwork_alpha_beta_args::m_in, "The number of edges being created in the last time interval.")
        .def_readwrite("m_out", &flockwork_alpha_beta_args::m_out, "The number of edges being deleted in the last time interval.")
        .def_readwrite("new_time", &flockwork_alpha_beta_args::new_time, "The bin edges of the new time bins");

    py::class_<edge_trajectories>(m, "edge_trajectories", R"pbdoc(
        Instead of getting lists of edges ordered in time, this description
        of a temporal network consists of a list of :class:`_tacoma.edge_trajectory_entry`. Each entry
        contains the edge and a list of time pairs denoting the times the edge exists.

        Optionally, dependent on the function which created this object, this object can
        contain edge similarities, which are defined as follows.
        Each edge :math:`i`, where
        :math:`i` is the edge's index in the trajectory list, is considered similar to edge :math:`j`
        if both are connected to the same node at the same time. Hence, if both :math:`i=(u,v)`
        and :math:`j=(u,w)` have
        node :math:`u` in common, their similarity is

        .. math::

            E_{ij} = \int\limits_0^{t_{\mathrm{max}}}dt\ A_{uv}(t)A_{uw}(t)
    )pbdoc")
        .def(py::init<>())
        .def_readwrite("trajectories", &edge_trajectories::trajectories,
                       R"pbdoc(Each entry of this list has properties `.edge` containing
                its nodes and `.time_pairs` containing the time intervals when
                the edge was active.)pbdoc")
        .def_readwrite("edge_similarities", &edge_trajectories::edge_similarities,
                       R"pbdoc(
                Each entry of this list is a triple `(i, j, w)`, where `i`
                is the `i`-th edge in `trajectories` (similarly for `j`) and `w`
                is their similarity.
                )pbdoc")
        .def_readwrite("int_to_node", &edge_trajectories::int_to_node,
                       R"pbdoc(A dictionary int -> string which keeps
                the original node names.)pbdoc")
        .def_readwrite("time_unit", &edge_trajectories::time_unit,
                       R"pbdoc(A string containing the unit of
                time for this network.)pbdoc")
        .def_readwrite("notes", &edge_trajectories::notes,
                       R"pbdoc(A string containing additional notes for this network.)pbdoc")
        .def_readwrite("N", &edge_trajectories::N,
                       R"pbdoc(Number of nodes)pbdoc")
        .def_readwrite("t0", &edge_trajectories::t0,
                       R"pbdoc(The initial time)pbdoc")
        .def_readwrite("tmax", &edge_trajectories::tmax,
                       R"pbdoc(The final time)pbdoc");

    py::class_<QS_SIS>(m, "QS_SIS", R"pbdoc(Base class for the simulation of an quasi-stationary SIS 
                                            compartmental infection model on a temporal network. Pass 
                                            this to :func:`tacoma.api.quasistationary_simulation` 
                                            to simulate and retrieve the simulation results.)pbdoc")
        .def(py::init<size_t, double, double, double, size_t, double, size_t, size_t, bool, size_t, bool>(),
             py::arg("N"),
             py::arg("t_simulation"),
             py::arg("infection_rate"),
             py::arg("recovery_rate"),
             py::arg("N_QS_samples"),
             py::arg("sampling_rate"),
             py::arg("number_of_initially_infected") = 1,
             py::arg("number_of_initially_vaccinated") = 0,
             py::arg("sample_network_state") = true,
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
                    N_QS_samples : int
                        Number of quasi-stationary configuration samples to be saved during the simulation.
                    sampling_rate : float
                        Rate with which to sample for the quasi-stationary configuration collection.
                    number_of_initially_infected : int, default = 1
                        Number of nodes which will be in the infected compartment at :math:`t = t_0`. Note that the default
                        value 1 is not an ideal initial value as fluctuations may lead to a quick end of the simulation
                        skewing the outcome. I generally recommend to use a number of the order of :math:`N/2`.
                    number_of_initially_vaccinated : int, default = 0
                        Number of nodes which will be in the recovered compartment at :math:`t = t_0`.
                    sample_network_state : bool, default = True
                        Do not only sample the node stati but also the current network structure.
                    seed : int, default = 0
                        Seed for RNG initialization. If this is 0, the seed will be initialized randomly.
                    verbose : bool, default = False
                        Be talkative.
                )pbdoc")
        .def_readwrite("configurations", &QS_SIS::QS_samples,
                R"pbdoc(Sampled collection. 
                )pbdoc"
            )
        .def_readwrite("last_active_time", &QS_SIS::last_active_time,
                R"pbdoc(The last time the model was active. 
                )pbdoc"
            )
        .def_readwrite("t_simulation", &QS_SIS::t_simulation,
                R"pbdoc(Time to simulate from :t0:.
                )pbdoc"
            )
        .def("ended_in_absorbing_state",
             &QS_SIS::simulation_ended,
             R"pbdoc(Return whether or not the simulation ended in an absorbing state.
             )pbdoc"
            )
        .def("get_random_configuration",
             &QS_SIS::get_random_configuration,
             R"pbdoc(Get a random configuration from the quasi-stationary
                     collection. Returns an `N`-length vector filled with 
                     the node stati and the Graph of the configuration.
             )pbdoc"
            )
        .def("set_initial_configuration",
             &QS_SIS::set_initial_configuration,
             R"pbdoc(Set a time and node statii.
             )pbdoc"
            )
        .def("get_infection_observables",
             &QS_SIS::get_infection_observables,
             R"pbdoc(Returns :math:`\left\langle I \right\rangle` and 
                     :math:`\left\langle I^2 \right\rangle` where the 
                     average is taken over the collection of quasi-stationary 
                     configurations.
             )pbdoc"
            );

    py::class_<MARKOV_SIS>(m, "MARKOV_SIS", "Base class for the markov integration of an SIS compartmental infection model on a temporal network. Pass this to :func:`tacoma.api.markov_SIS` to simulate and retrieve the simulation results.")
        .def(py::init<size_t, double, double, double, double, size_t, size_t, double, size_t, bool>(),
             py::arg("N"),
             py::arg("t_simulation"),
             py::arg("infection_rate"),
             py::arg("recovery_rate"),
             py::arg("minimum_I"),
             py::arg("number_of_initially_infected") = 1,
             py::arg("number_of_initially_vaccinated") = 0,
             py::arg("sampling_dt") = 0.0,
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
                    minimum_I : float
                        If the total sum of infection probability over all nodes is below this threshold, the disease
                        is considered to be died out.
                    number_of_initially_infected : int, default = 1
                        Number of nodes which will be in the infected compartment at :math:`t = t_0`. Note that the default
                        value 1 is not an ideal initial value as fluctuations may lead to a quick end of the simulation
                        skewing the outcome. I generally recommend to use a number of the order of :math:`N/2`.
                    number_of_initially_vaccinated : int, default = 0
                        Number of nodes which will be in the recovered compartment at :math:`t = t_0`.
                    sampling_dt : float, default = 0.0
                        If this is ``>0.0``, save observables roughly every sampling_dt instead of on every change.
                    seed : int, default = 0
                        Seed for RNG initialization. If this is 0, the seed will be initialized randomly.
                    verbose : bool, default = False
                        Be talkative.
                )pbdoc")
        .def_readwrite("time", &MARKOV_SIS::time, "A list containing the time points at which one or more of the observables changed.")
        .def_readwrite("R0", &MARKOV_SIS::R0, R"pbdoc(
                   A list containing the basic reproduction number defined as :math:`R_0(t) = \eta\left\langle k \right\rangle(t) / \rho`
                   where :math:`\eta` is the infection rate per link and :math:`\rho` is the recovery rate per node.
                   )pbdoc")
        .def_readwrite("I", &MARKOV_SIS::I, "A list containing accumulated probability to be infected at time :math:`t`.")
        .def_readwrite("t_simulation", &MARKOV_SIS::t_simulation, "Absolute run time of the simulation.");

    py::class_<SIS>(m, "SIS", "Base class for the simulation of an SIS compartmental infection model on a temporal network. Pass this to :func:`tacoma.api.gillespie_SIS` to simulate and retrieve the simulation results.")
        .def(py::init<size_t, double, double, double, size_t, size_t, bool, bool, double, size_t, bool>(),
             py::arg("N"),
             py::arg("t_simulation"),
             py::arg("infection_rate"),
             py::arg("recovery_rate"),
             py::arg("number_of_initially_infected") = 1,
             py::arg("number_of_initially_vaccinated") = 0,
             py::arg("prevent_disease_extinction") = false,
             py::arg("save_infected_nodes") = false,
             py::arg("sampling_dt") = 0.0,
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
                    prevent_disease_extinction : bool, default = False
                        If this is `True`, the recovery of the last infected node will always be prevented.
                    save_infected_nodes : bool, default = False
                        If this is `True`, the recovery of the last infected node will always be prevented.
                    sampling_dt : float, default = 0.0
                        If this is ``>0.0``, save observables roughly every sampling_dt instead of on every change.
                    seed : int, default = 0
                        Seed for RNG initialization. If this is 0, the seed will be initialized randomly.
                    verbose : bool, default = False
                        Be talkative.
                )pbdoc")
        .def("set_node_status",&SIS::set_node_status,
             R"pbdoc(Set the state of the epidemic as given by a list :obj:`list` of :obj:`int`.
                     Needs to have the same lengths as the number of nodes :math:`N`.)pbdoc")
        .def("get_node_status",&SIS::get_node_status,
             R"pbdoc(Retrieve the current state of the epidemic as given by a list :obj:`list` of :obj:`int`. Will have the same lengths as the number of nodes :math:`N`, each entry `i` containing the status of node `i`.)pbdoc")
        .def_readwrite("time", &SIS::time, "A list containing the time points at which one or more of the observables changed.")
        .def_readwrite("R0", &SIS::R0, R"pbdoc(
                   A list containing the basic reproduction number defined as :math:`R_0(t) = \eta\left\langle k \right\rangle(t) / \rho`
                   where :math:`\eta` is the infection rate per link and :math:`\rho` is the recovery rate per node.
                   )pbdoc")
        .def_readwrite("SI", &SIS::SI, "A list containing the number of :math:`SI`-links at time :math:`t`.")
        .def_readwrite("I", &SIS::I, "A list containing the number of infected at time :math:`t`.")
        .def_readwrite("infected_nodes", &SIS::saved_infected_nodes, "A list of lists of ints. Each list of ints corresponds to the node integers which were infected at the corresonding time :math:`t`.")
        .def_readwrite("t_simulation", &SIS::t_simulation, "Absolute run time of the simulation.");

    py::class_<cluster_size_SIS>(m, "cluster_size_SIS", R"pbdoc(Base class for the simulation of an SIS compartmental infection model on a temporal network. Pass this to :func:`tacoma.api.gillespie_epidemics` to simulate and retrieve the simulation results. This simulation runs until all initially infected nodes recovered at least once, which is useful to find the epidemic threshold by measuring the impact of the seed(s). If you want to sample the standard observables, set ``sampling_dt>=0.0``.)pbdoc")
        .def(py::init<size_t, double, double, double, size_t, size_t, double, size_t, bool>(),
             py::arg("N"),
             py::arg("t_simulation"),
             py::arg("infection_rate"),
             py::arg("recovery_rate"),
             py::arg("number_of_initially_infected") = 1,
             py::arg("number_of_initially_vaccinated") = 0,
             py::arg("sampling_dt") = -1.0,
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
                    sampling_dt : float, default = -1.0
                        If it is negative, do not save any observables but the life time of the process.
                        If this is ``>0.0``, save observables roughly every sampling_dt instead of on every change.
                        If it is 0.0, sample at every change.
                    seed : int, default = 0
                        Seed for RNG initialization. If this is 0, the seed will be initialized randomly.
                    verbose : bool, default = False
                        Be talkative.
                )pbdoc")
        .def_readwrite("time", &cluster_size_SIS::time, "A list containing the time points at which one or more of the observables changed.")
        .def_readwrite("lifetime", &cluster_size_SIS::lifetime, "The lifetime of the process.")
        .def_readwrite("coverage", &cluster_size_SIS::coverage, "The total number of nodes which have been infected at least once when the simulation ends.")
        .def_readwrite("cluster_size", &cluster_size_SIS::cluster_size, "The number of nodes which are infected at the time when the simulation ends.")
        .def_readwrite("number_of_events", &cluster_size_SIS::number_of_events, "The number of events happened during the process.")
        .def_readwrite("R0", &cluster_size_SIS::R0, R"pbdoc(
                   A list containing the basic reproduction number defined as :math:`R_0(t) = \eta\left\langle k \right\rangle(t) / \rho`
                   where :math:`\eta` is the infection rate per link and :math:`\rho` is the recovery rate per node.
                   )pbdoc")
        .def_readwrite("SI", &cluster_size_SIS::SI, "A list containing the number of :math:`SI`-links at time :math:`t`.")
        .def_readwrite("I", &cluster_size_SIS::I, "A list containing the number of infected at time :math:`t`.")
        .def_readwrite("t_simulation", &cluster_size_SIS::t_simulation, "Absolute run time of the simulation.");

    py::class_<coverage_SIS>(m, "coverage_SIS", R"pbdoc(Base class for the simulation of an SIS compartmental infection model on a temporal network. Pass this to :func:`tacoma.api.gillespie_SIS` to simulate and retrieve the simulation results. This simulation runs until a certain amount of nodes have been infected at least once, which is useful to fifind the epidemic threshold and thus, only the lifetime of the process is measured. If you want to sample the standard observables, set ``sampling_dt>=0.0``.)pbdoc")
        .def(py::init<size_t, double, double, double, size_t, size_t, double, double, size_t, bool>(),
             py::arg("N"),
             py::arg("t_simulation"),
             py::arg("infection_rate"),
             py::arg("recovery_rate"),
             py::arg("number_of_initially_infected") = 1,
             py::arg("number_of_initially_vaccinated") = 0,
             py::arg("critical_coverage") = 0.5,
             py::arg("sampling_dt") = -1.0,
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
                    critical_coverage : float, default = 0.5
                        The simulation ends if the total number of nodes which were infected at least once exceeds
                        this ratio.
                    sampling_dt : float, default = -1.0
                        If it is negative, do not save any observables but the life time of the process.
                        If this is ``>0.0``, save observables roughly every sampling_dt instead of on every change.
                        If it is 0.0, sample at every change.
                    seed : int, default = 0
                        Seed for RNG initialization. If this is 0, the seed will be initialized randomly.
                    verbose : bool, default = False
                        Be talkative.
                )pbdoc")
        .def_readwrite("time", &coverage_SIS::time, "A list containing the time points at which one or more of the observables changed.")
        .def_readwrite("lifetime", &coverage_SIS::lifetime, "The lifetime of the process.")
        .def_readwrite("number_of_events", &coverage_SIS::number_of_events, "The number of events happened during the process.")
        .def_readwrite("R0", &coverage_SIS::R0, R"pbdoc(
                   A list containing the basic reproduction number defined as :math:`R_0(t) = \eta\left\langle k \right\rangle(t) / \rho`
                   where :math:`\eta` is the infection rate per link and :math:`\rho` is the recovery rate per node.
                   )pbdoc")
        .def_readwrite("SI", &coverage_SIS::SI, "A list containing the number of :math:`SI`-links at time :math:`t`.")
        .def_readwrite("I", &coverage_SIS::I, "A list containing the number of infected at time :math:`t`.")
        .def_readwrite("t_simulation", &coverage_SIS::t_simulation, "Absolute run time of the simulation.");

    py::class_<eSIS>(m, "eSIS", R"pbdoc(Base class for the simulation of an :math:`\varepsilon`-SIS compartmental infection model on a temporal network. Pass this to :func:`tacoma.api.gillespie_epidemics` to simulate and retrieve the simulation results.)pbdoc")
        .def(py::init<size_t, double, double, double, double, size_t, size_t, bool, double, size_t, bool>(),
             py::arg("N"),
             py::arg("t_simulation"),
             py::arg("infection_rate"),
             py::arg("recovery_rate"),
             py::arg("self_infection_rate"),
             py::arg("number_of_initially_infected") = 1,
             py::arg("number_of_initially_vaccinated") = 0,
             py::arg("prevent_disease_extinction") = false,
             py::arg("sampling_dt") = 0.0,
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
                    self_infection_rate : float
                        Infection rate per susecptible (expected number of reaction events :math:`S\rightarrow I`
                        for a single susceptible node per dimension of time).
                    number_of_initially_infected : int, default = 1
                        Number of nodes which will be in the infected compartment at :math:`t = t_0`. Note that the default
                        value 1 is not an ideal initial value as fluctuations may lead to a quick end of the simulation
                        skewing the outcome. I generally recommend to use a number of the order of :math:`N/2`.
                    number_of_initially_vaccinated : int, default = 0
                        Number of nodes which will be in the recovered compartment at :math:`t = t_0`.
                    prevent_disease_extinction : bool, default = False
                        If this is `True`, the recovery of the last infected node will always be prevented.
                    sampling_dt : float, default = 0.0
                        If this is ``>0.0``, save observables roughly every sampling_dt instead of on every change.
                    seed : int, default = 0
                        Seed for RNG initialization. If this is 0, the seed will be initialized randomly.
                    verbose : bool, default = False
                        Be talkative.
                )pbdoc")
        .def_readwrite("time", &eSIS::time, "A list containing the time points at which one or more of the observables changed.")
        .def_readwrite("R0", &eSIS::R0, R"pbdoc(
                   A list containing the basic reproduction number defined as :math:`R_0(t) = \eta\left\langle k \right\rangle(t) / \rho`
                   where :math:`\eta` is the infection rate per link and :math:`\rho` is the recovery rate per node.
                   )pbdoc")
        .def_readwrite("SI", &eSIS::SI, "A list containing the number of :math:`SI`-links at time :math:`t`.")
        .def_readwrite("I", &eSIS::I, "A list containing the number of infected at time :math:`t`.")
        .def_readwrite("t_simulation", &eSIS::t_simulation, "Absolute run time of the simulation.");

    py::class_<node_based_SIS>(m, "node_based_SIS", R"pbdoc(Base class for the simulation of an SIS compartmental infection model on a temporal network using an SI-Graph for keeping track of SI-edges (meaning this is a node-based algorithm). Pass this to :func:`tacoma.api.gillespie_node_based_SIS` to simulate and retrieve the simulation results.)pbdoc")
        .def(py::init<size_t, double, double, double, size_t, size_t, bool, double, size_t, bool>(),
             py::arg("N"),
             py::arg("t_simulation"),
             py::arg("infection_rate"),
             py::arg("recovery_rate"),
             py::arg("number_of_initially_infected") = 1,
             py::arg("number_of_initially_vaccinated") = 0,
             py::arg("prevent_disease_extinction") = false,
             py::arg("sampling_dt") = 0.0,
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
                    sampling_dt : float, default = 0.0
                        If this is ``>0.0``, save observables roughly every sampling_dt instead of on every change.
                    seed : int, default = 0
                        Seed for RNG initialization. If this is 0, the seed will be initialized randomly.
                    verbose : bool, default = False
                        Be talkative.
                )pbdoc")
        .def_readwrite("time", &node_based_SIS::time, "A list containing the time points at which one or more of the observables changed.")
        .def_readwrite("R0", &node_based_SIS::R0, R"pbdoc(
                   A list containing the basic reproduction number defined as :math:`R_0(t) = \eta\left\langle k \right\rangle(t) / \rho`
                   where :math:`\eta` is the infection rate per link and :math:`\rho` is the recovery rate per node.
                   )pbdoc")
        .def_readwrite("SI", &node_based_SIS::SI, "A list containing the number of :math:`SI`-links at time :math:`t`.")
        .def_readwrite("I", &node_based_SIS::I, "A list containing the number of infected at time :math:`t`.")
        .def_readwrite("t_simulation", &node_based_SIS::t_simulation, "Absolute run time of the simulation.");


    py::class_<SI>(m, "SI", "Base class for the simulation of an SI compartmental infection model on a temporal network. Pass this to :func:`tacoma.api.gillespie_SI` to simulate and retrieve the simulation results. Simulation stops when ``t_simulation`` is reached or if no infected is left.")
        .def(py::init<size_t, double, double, size_t, size_t, double, size_t, bool, bool>(),
             py::arg("N"),
             py::arg("t_simulation"),
             py::arg("infection_rate"),
             py::arg("number_of_initially_infected") = 1,
             py::arg("number_of_initially_vaccinated") = 0,
             py::arg("sampling_dt") = 0.0,
             py::arg("seed") = 0,
             py::arg("save_infection_events") = false,
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
                    sampling_dt : float, default = 0.0
                        If this is ``>0.0``, save observables roughly every sampling_dt instead of on every change.
                    seed : int, default = 0
                        Seed for RNG initialization. If this is 0, the seed will be initialized randomly.
                    save_infection_events: bool, default = False
                        If true, the edge along which each infection event occurs is saved in the variable `infection_events`.
                    verbose : bool, default = False
                        Be talkative.
                )pbdoc")
        .def_readwrite("time", &SI::time, "A list containing the time points at which one or more of the observables changed.")
        .def_readwrite("SI", &SI::_SI, "A list containing the number of :math:`SI`-links at time :math:`t`.")
        .def_readwrite("I", &SI::I, "A list containing the number of infected at time :math:`t`.")
        .def_readwrite("infection_events", &SI::infection_events, "A list containing the edges along which each infection event took place, in the form (infection_source, susceptible).")
        .def_readwrite("t_simulation", &SI::t_simulation, "Absolute run time of the simulation.");

    py::class_<SIR>(m, "SIR", "Base class for the simulation of an SIR compartmental infection model on a temporal network. Pass this to :func:`tacoma.api.gillespie_SIR` to simulate and retrieve the simulation results.")
        .def(py::init<size_t, double, double, double, size_t, size_t, double, size_t, bool>(),
             py::arg("N"),
             py::arg("t_simulation"),
             py::arg("infection_rate"),
             py::arg("recovery_rate"),
             py::arg("number_of_initially_infected") = 1,
             py::arg("number_of_initially_vaccinated") = 0,
             py::arg("sampling_dt") = 0.0,
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
                    sampling_dt : float, default = 0.0
                        If this is ``>0.0``, save observables roughly every sampling_dt instead of on every change.
                    seed : int, default = 0
                        Seed for RNG initialization. If this is 0, the seed will be initialized randomly.
                    verbose : bool, default = False
                        Be talkative.
                )pbdoc")
        .def_readwrite("time", &SIR::time, "A list containing the time points at which one or more of the observables changed.")
        .def_readwrite("R0", &SIR::R0, R"pbdoc(
                   A list containing the basic reproduction number defined as :math:`R_0(t) = \eta\left\langle k \right\rangle(t) / \rho`
                   where :math:`\eta` is the infection rate per link and :math:`\rho` is the recovery rate per node.
                   )pbdoc")
        .def_readwrite("SI", &SIR::SI, "A list containing the number of :math:`SI`-links at time :math:`t`.")
        .def_readwrite("I", &SIR::I, "A list containing the number of infected at time :math:`t`.")
        .def_readwrite("R", &SIR::R, "A list containing the number of recovered at time :math:`t`.")
        .def_readwrite("t_simulation", &SIR::t_simulation, "Absolute run time of the simulation.");

    py::class_<SIRS>(m, "SIRS", "Base class for the simulation of an SIRS compartmental infection model on a temporal network. Pass this to :func:`tacoma.api.gillespie_SIRS` to simulate and retrieve the simulation results.")
        .def(py::init<size_t, double, double, double, double, size_t, size_t, double, size_t, bool>(),
             py::arg("N"),
             py::arg("t_simulation"),
             py::arg("infection_rate"),
             py::arg("recovery_rate"),
             py::arg("waning_immunity_rate"),
             py::arg("number_of_initially_infected") = 1,
             py::arg("number_of_initially_vaccinated") = 0,
             py::arg("sampling_dt") = 0.0,
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
                        Number of nodes which will be in the vaccinated compartment at :math:`t = t_0`.
                    sampling_dt : float, default = 0.0
                        If this is ``>0.0``, save observables roughly every sampling_dt instead of on every change.
                    seed : int, default = 0
                        Seed for RNG initialization. If this is 0, the seed will be initialized randomly.
                    verbose : bool, default = False
                        Be talkative.
                )pbdoc")
        .def_readwrite("time", &SIRS::time, "A list containing the time points at which one or more of the observables changed.")
        .def_readwrite("R0", &SIRS::R0, R"pbdoc(
                    A list containing the basic reproduction number defined as :math:`R_0(t) = \eta\left\langle k \right\rangle(t) / \rho`
                    where :math:`\eta` is the infection rate per link and :math:`\rho` is the recovery rate per node.
            )pbdoc")
        .def_readwrite("SI", &SIRS::SI, "A list containing the number of :math:`SI`-links at time :math:`t`.")
        .def_readwrite("I", &SIRS::I, "A list containing the number of infected at time :math:`t`.")
        .def_readwrite("R", &SIRS::R, "A list containing the number of recovered at time :math:`t`.")
        .def_readwrite("t_simulation", &SIRS::t_simulation, "Absolute run time of the simulation.");

    py::class_<EdgeActivityModel>(m, "EdgeActivityModel",
            R"pbdoc(
                Base class for the simulation of a simple edge activity model. Pass this to :func:`tacoma.api.gillespie_epidemics` or
                :func:`tacoma.api.markov_epidemics`.
            )pbdoc")
        .def(py::init<size_t, double, double, double, bool, bool, size_t, bool>(),
             py::arg("N"),
             py::arg("rho"),
             py::arg("omega"),
             py::arg("t0") = 0.0,
             py::arg("use_rejection_sampling_of_non_neighbor") = true,
             py::arg("save_temporal_network") = false,
             py::arg("seed") = 0,
             py::arg("verbose") = false,
             R"pbdoc(
                    Parameters
                    ----------
                    N : int
                        Number of nodes in the temporal network.
                    rho : float
                        Demanded network density.
                    omega : float
                        rate with which edges are switched on and off, respectively,
                        :math:`\omega^{-1}=(\omega^-)^{-1} + (\omega^+)^{-1}`.
                    t0 : float, default = 0.0
                        initial time
                    use_rejection_sampling_of_non_neighbor : bool, default: True
                        If this is `True`, the edges to be turned on are sampled 
                        by drawing random edges until one is found which is turned
                        off. If `False`, there's a more sophisticated but probably
                        slower method (use this option for dense networks).
                    save_temporal_network : bool, default: False
                        If this is `True`, the changes are saved in an instance of 
                        :func:`_tacoma.edge_changes` (in the attribute `edge_changes`.
                    seed : int, default = 0
                        Seed for RNG initialization. If this is 0, the seed will be initialized randomly.
                        However, the generator will be rewritten 
                        in :func:`tacoma.api.gillespie_SIS_EdgeActivityModel` anyway.
                    verbose : bool, default = False
                        Be talkative.
                )pbdoc")
        .def("set_initial_configuration",&EdgeActivityModel::set_initial_configuration,
             R"pbdoc(Reset the state of the network to a certain graph (:obj:`list` of :obj:`set` of :obj:`int`))pbdoc")
        .def("set_initial_edgelist",&EdgeActivityModel::set_initial_edgelist,
             R"pbdoc(Reset the state of the network to a certain edgelist (:obj:`list` of :obj:`tuple` of :obj:`int`))pbdoc")
        .def("get_current_edgelist",&EdgeActivityModel::get_current_edgelist,
             R"pbdoc(Get an edge list of the current network state.)pbdoc")
        .def_readwrite("edge_changes", &EdgeActivityModel::edg_chg, 
                    R"pbdoc(An instance of :class:`_tacoma.edge_changes` with the saved temporal network (only if
                    `save_temporal_network` is `True`).)pbdoc")
        .def_readwrite("N", &EdgeActivityModel::N, 
                    R"pbdoc(Number of nodes.)pbdoc");

    py::class_<FlockworkPModel>(m, "FlockworkPModel", 
            R"pbdoc(
                Base class for the simulation of a simple Flockwork-P model. Pass this to :func:`tacoma.api.gillespie_epidemics` or
                :func:`tacoma.api.markov_epidemics`.
            )pbdoc")
        .def(py::init< vector < pair < size_t, size_t > > , size_t, double, double, double, bool, size_t, bool>(),
             py::arg("E"),
             py::arg("N"),
             py::arg("gamma"),
             py::arg("P"),
             py::arg("t0") = 0.0,
             py::arg("save_temporal_network") = false,
             py::arg("seed") = 0,
             py::arg("verbose") = false,
             R"pbdoc(
                    Parameters
                    ----------
                    E : list of pair of int
                        Initial edge list.
                    N : int
                        Number of nodes in the temporal network.
                    gamma : float
                        The probability per unit time per node that any event happens.
                    P : float
                        The probability to reconnect when an event happened.
                    t0 : float, default = 0.0
                        initial time
                    save_temporal_network : bool, default: False
                        If this is `True`, the changes are saved in an instance of 
                        :func:`_tacoma.edge_changes` (in the attribute `edge_changes`.
                    seed : int, default = 0
                        Seed for RNG initialization. If this is 0, the seed will be initialized randomly.
                        However, the generator will be rewritten 
                        in :func:`tacoma.api.gillespie_SIS_EdgeActivityModel` anyway.
                    verbose : bool, default = False
                        Be talkative.
                )pbdoc")
        .def("set_initial_configuration",&FlockworkPModel::set_initial_configuration,
             R"pbdoc(Reset the state of the network to a certain graph (:obj:`list` of :obj:`set` of :obj:`int`))pbdoc")
        .def("get_current_edgelist",&FlockworkPModel::get_current_edgelist,
             R"pbdoc(Get an edge list of the current network state.)pbdoc")
        .def("simulate",&FlockworkPModel::simulate,
             py::arg("t_run_total"), 
             py::arg("reset") = true,
             py::arg("save_temporal_network") = true,
             R"pbdoc(Simulate a Flockwork model until ``t_run_total``.
                     Obtain the result via ``FlockworkPModel.edge_changes``.)pbdoc")
        .def_readwrite("edge_changes", &FlockworkPModel::edg_chg, 
                    R"pbdoc(An instance of :class:`_tacoma.edge_changes` with the saved temporal network (only if
                    `save_temporal_network` is `True`).)pbdoc")
        .def_readwrite("N", &FlockworkPModel::N, 
                    R"pbdoc(Number of nodes.)pbdoc");

}
