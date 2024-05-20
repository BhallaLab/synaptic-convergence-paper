#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : plot_figure3.py
# Author            : Bhanu Priya S <bhanups@ncbs.res.in>
# Date              : 09.11.2023
# Last Modified Date: 09.11.2023
# Last Modified By  : Bhanu Priya S <bhanups@ncbs.res.in>

"""Plot Figure-3
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.image as mpimg
import matplotlib.patches as mpatches
import seaborn as sns
import pandas as pd
from scipy.stats import poisson
from scipy.linalg import expm

sys.path.append("../../codes/sequences/")
import sequence_params as prm

sns.set(style="ticks")
sns.set_context("paper")

DATA_PATH = "../../codes/sequences/data/"

networks = [
    prm.hippo_chem,
    prm.hippo_cicr,
    prm.hippo_elec,
    prm.cortex_chem,
    prm.cortex_cicr,
    prm.cortex_elec,
]
colors = sns.color_palette("muted", n_colors=len(networks))

MIN_M = prm.MIN_M
MAX_M = prm.MAX_M
T_POP_POST = 4000
NUM_RUNS = 100


def read_data(net, prefix, trial_numbers):
    """Read sequence counts data across trials

    Parameters
    ----------
    net : dict
        Dictionary of network parameters
    prefix : str
        Type of data file to be read
    trial_numbers : np.ndarray
        Array of trial numbers to read

    Returns
    -------
    data_matrix : np.ndarray
        3D array of sequence counts of shape num_trials x (t_pop_post x num_runs) x num_sequence_sizes
    """
    data_matrix = np.zeros(
        (len(trial_numbers), T_POP_POST * NUM_RUNS, MAX_M - MIN_M + 1), dtype=np.int32
    )
    for t, t_num in enumerate(trial_numbers):
        for r in range(NUM_RUNS):
            filename = DATA_PATH + "%s-sequences-%s-run-%03d-trial-%03d.csv" % (
                prefix,
                net["name"],
                r + 1,
                t_num,
            )
            data = np.asarray(pd.read_csv(filename, header=None))
            data_matrix[t, r * T_POP_POST : (r + 1) * T_POP_POST, :] = data

    return data_matrix


def plot_sequence_types(ax):
    """Display schematic of gap-fill and noise sequences

    Parameters
    ----------
    ax : axis object
    """
    ax.imshow(mpimg.imread("sequence_types.png"))
    ax.text(
        -0.04,
        0.13,
        "Gap-fill\nsequence",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.text(
        -0.08,
        0.72,
        "Perfect\nstimulus-driven\nsequence",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.text(
        0.95,
        0.19,
        "Noise\nsequence",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.axis("off")


def plot_active_vs_inactive_sequences(ax):
    """Distinguishing active sequences from inactive sequences

    Parameters
    ----------
    ax : axis object
    """
    ax.imshow(mpimg.imread("active_vs_inactive_sequences.png"))
    ax.text(
        -0.04,
        0.13,
        "Active\nsequence",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.text(
        0.95,
        0.19,
        "Inactive\nsequence",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.axis("off")


def plot_sequences_venn_diagram(ax, bx):
    """Distinguishing the venn diagram of sequences

    Parameters
    ----------
    ax : axis object
    bx : axis object
    """
    ax.imshow(mpimg.imread("sequences_venn_diagram.png"))

    ax.plot([0.17, 0.11], [0.88, 0.91], color="k", transform=ax.transAxes)
    ax.text(
        0.08,
        0.92,
        "U",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.plot([0.58, 0.79], [0.55, 0.73], color="k", transform=ax.transAxes)
    ax.text(
        0.85,
        0.75,
        "PSASD",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.plot([0.65, 0.79], [0.39, 0.28], color="k", transform=ax.transAxes)
    ax.text(
        0.85,
        0.25,
        "PSCSD",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.axis("off")
    legend_handles = []
    neuron_types_colors = [
        "#000000FF",
        "#AA4499FF",
        "#8877DEFF",
    ]
    neuron_types = [
        "$U$ - All postsynaptic neurons",
        "$PSCSD$ - Neurons receiving perfect sequential connections from stimulus-driven neurons",
        "$PSASD$ - Neurons receiving perfect sequential activity from stimulus-driven neurons",
    ]
    for i, nt in enumerate(neuron_types):
        patch = mpatches.Patch(
            edgecolor=neuron_types_colors[i], facecolor="none", label=nt
        )
        legend_handles.append(patch)
    bx.legend(handles=legend_handles, loc="center", frameon=False)
    bx.axis("off")


def plot_perfect_stimulus_driven_conn_sequences(ax):
    """Plot probability of true sequences based on connectivity

    Parameters
    ----------
    ax : axis object
    """
    for i, net in enumerate(networks):
        stim_array = np.genfromtxt(
            DATA_PATH + f"stimulus-{net['name']}-run-001.csv", delimiter=","
        )
        stim_trial_numbers = np.where(stim_array < 0)[0] + 1
        perfect_stimulus_driven_sequence_counts = read_data(
            net, "perfect-stimulus-driven", stim_trial_numbers
        )
        p_obs = np.mean(perfect_stimulus_driven_sequence_counts > 0, axis=(0, 1))

        p = net["connection_prob"]
        N = net["ensemble_size"]
        L = net["dendrite_length"]
        M = np.arange(MIN_M, MAX_M + 1)
        delta = net["delta"]

        E_seq_exp = p * N * np.power(p * N * delta / L, M - 1)
        p_seq_exp = 1 - np.exp(-E_seq_exp)
        ax.plot(M, p_seq_exp, ":x", c=colors[i], label=net["name"])
        ax.plot(M, p_obs, "-o", c=colors[i], label=net["name"])

    sns.despine()
    ax.set_yscale("log", nonpositive="mask")
    ax.set_ylim([1e-5, 1e0])
    ax.set_ylabel("p_cposs")
    ax.set_xticks([3, 5, 7, 9])
    ax.text(
        0.5,
        0.9,
        "Connectivity based",
        horizontalalignment="center",
        verticalalignment="center",
        fontsize="large",
        transform=ax.transAxes,
    )


def plot_perfect_stimulus_driven_active_sequences(ax):
    """Plot probability of true sequences based on activity

    Parameters
    ----------
    ax : axis object
    """
    for i, net in enumerate(networks):
        stim_array = np.genfromtxt(
            DATA_PATH + f"stimulus-{net['name']}-run-001.csv", delimiter=","
        )
        stim_trial_numbers = np.where(stim_array > 0)[0] + 1
        perfect_stimulus_driven_sequence_counts = read_data(
            net, "perfect-stimulus-driven", stim_trial_numbers
        )
        p_obs = np.mean(perfect_stimulus_driven_sequence_counts > 0, axis=(0, 1))

        p = net["connection_prob"]
        p_e = net["ensemble_participation_prob"]
        N = net["ensemble_size"]
        L = net["dendrite_length"]
        M = np.arange(MIN_M, MAX_M + 1)
        delta = net["delta"]

        E_seq_exp = p * p_e * N * np.power(p * p_e * N * delta / L, M - 1)
        p_seq_exp = 1 - np.exp(-E_seq_exp)
        ax.plot(M, p_seq_exp, ":x", c=colors[i], label=net["name"])
        ax.plot(M, p_obs, "-o", c=colors[i], label=net["name"])

    sns.despine()
    ax.set_yscale("log", nonpositive="mask")
    ax.set_ylim([1e-5, 1e0])
    ax.set_ylabel("p_aposs")
    ax.set_xticks([3, 5, 7, 9])
    ax.set_xlabel("Sequence length M")
    ax.text(
        0.5,
        0.9,
        "Activity based",
        horizontalalignment="center",
        verticalalignment="center",
        fontsize="large",
        transform=ax.transAxes,
    )


def plot_noise_sequences(ax):
    """Plot probability of noise sequences based on activity

    Parameters
    ----------
    ax : axis object
    """
    for i, net in enumerate(networks):
        stim_array = np.genfromtxt(
            DATA_PATH + f"stimulus-{net['name']}-run-001.csv", delimiter=","
        )
        stim_trial_numbers = np.where(stim_array > 0)[0] + 1
        noise_sequence_counts = read_data(net, "noise", stim_trial_numbers)
        p_obs = np.mean(noise_sequence_counts > 0, axis=(0, 1))

        p = net["connection_prob"]
        N = net["ensemble_size"]
        L = net["dendrite_length"]
        M = np.arange(MIN_M, MAX_M + 1)
        delta = net["delta"]
        t_pop_pre = net["population_size"]
        R = net["background_rate"]
        D = net["sequence_time_step"]
        noise_prob = 1 - np.exp(-R * D)
        num_synapses = int(p * t_pop_pre)

        E_seq_exp = (
            noise_prob
            * (num_synapses - p * N)
            * np.power(noise_prob * (num_synapses - p * N) * delta / L, M - 1)
        )
        p_seq_exp = 1 - np.exp(-E_seq_exp)
        ax.plot(M, p_seq_exp, ":x", c=colors[i], label=net["name"])
        ax.plot(M, p_obs, "-o", c=colors[i], label=net["name"])

    sns.despine()
    ax.set_yscale("log", nonpositive="mask")
    ax.set_ylim([1e-5, 1e0])
    ax.set_ylabel("p_noise_seq")
    ax.set_xticks([3, 5, 7, 9])
    ax.set_yticklabels([])


def plot_any_sequences(ax):
    """Plot probability of any sequence based on activity

    Parameters
    ----------
    ax : axis object
    """
    for i, net in enumerate(networks):
        stim_array = np.genfromtxt(
            DATA_PATH + f"stimulus-{net['name']}-run-001.csv", delimiter=","
        )
        stim_trial_numbers = np.where(stim_array > 0)[0] + 1
        all_sequence_counts = read_data(net, "all", stim_trial_numbers)
        p_obs = np.mean(all_sequence_counts > 0, axis=(0, 1))

        p = net["connection_prob"]
        p_e = net["ensemble_participation_prob"]
        N = net["ensemble_size"]
        L = net["dendrite_length"]
        M = np.arange(MIN_M, MAX_M + 1)
        delta = net["delta"]
        t_pop_pre = net["population_size"]
        R = net["background_rate"]
        D = net["sequence_time_step"]
        noise_prob = 1 - np.exp(-R * D)
        num_synapses = int(p * t_pop_pre)

        E_seq_exp = (p * p_e * N + noise_prob * (num_synapses - p * N)) * np.power(
            (p * p_e * N + noise_prob * (num_synapses - p * N)) * delta / L, M - 1
        )
        p_seq_exp = 1 - np.exp(-E_seq_exp)
        ax.plot(M, p_seq_exp, ":x", c=colors[i], label=net["name"])
        ax.plot(M, p_obs, "-o", c=colors[i], label=net["name"])

    sns.despine()
    ax.set_yscale("log", nonpositive="mask")
    ax.set_ylim([1e-5, 1e0])
    ax.set_xlabel("Sequence length M")
    ax.set_ylabel("p_any_seq")
    ax.set_xticks([3, 5, 7, 9])
    ax.set_yticklabels([])


def plot_false_positive_sequences(ax):
    """Plot probability of false positive sequences

    Parameters
    ----------
    ax : axis object
    """
    for i, net in enumerate(networks):
        stim_array = np.genfromtxt(
            DATA_PATH + f"stimulus-{net['name']}-run-001.csv", delimiter=","
        )
        stim_trial_numbers = np.where(stim_array > 0)[0] + 1
        perfect_stimulus_driven_sequence_counts = read_data(
            net, "perfect-stimulus-driven", stim_trial_numbers
        )
        all_sequence_counts = read_data(net, "all", stim_trial_numbers)
        E_perfect_stimulus_driven_obs = np.mean(
            perfect_stimulus_driven_sequence_counts, axis=(0, 1)
        )
        E_all_obs = np.mean(all_sequence_counts, axis=(0, 1))
        E_fp_obs = E_all_obs - E_perfect_stimulus_driven_obs
        p_fp_obs = 1 - np.exp(-E_fp_obs)

        p = net["connection_prob"]
        p_e = net["ensemble_participation_prob"]
        N = net["ensemble_size"]
        L = net["dendrite_length"]
        M = np.arange(MIN_M, MAX_M + 1)
        delta = net["delta"]
        t_pop_pre = net["population_size"]
        R = net["background_rate"]
        D = net["sequence_time_step"]
        noise_prob = 1 - np.exp(-R * D)
        num_synapses = int(p * t_pop_pre)

        E_perfect_stimulus_driven_seq_exp = (
            p * p_e * N * np.power(p * p_e * N * delta / L, M - 1)
        )
        E_all_seq_exp = (p * p_e * N + noise_prob * (num_synapses - p * N)) * np.power(
            (p * p_e * N + noise_prob * (num_synapses - p * N)) * delta / L, M - 1
        )
        E_fp_seq_exp = E_all_seq_exp - E_perfect_stimulus_driven_seq_exp
        p_fp_seq_exp = 1 - np.exp(-E_fp_seq_exp)
        ax.plot(M, p_fp_seq_exp, ":x", c=colors[i], label=net["name"])
        ax.plot(M, p_fp_obs, "-o", c=colors[i], label=net["name"])

    sns.despine()
    ax.set_yscale("log", nonpositive="mask")
    ax.set_ylim([1e-5, 1e0])
    ax.set_xlabel("Sequence length M")
    ax.set_ylabel("p_false_positive")
    ax.set_xticks([3, 5, 7, 9])
    ax.set_yticklabels([])


def plot_gap_fill_sequences(ax):
    """Plot probability of gap fill sequences

    Parameters
    ----------
    ax : axis object
    """
    for i, net in enumerate(networks):
        stim_array = np.genfromtxt(
            DATA_PATH + f"stimulus-{net['name']}-run-001.csv", delimiter=","
        )
        stim_trial_numbers = np.where(stim_array > 0)[0] + 1
        perfect_stimulus_driven_sequence_counts = read_data(
            net, "perfect-stimulus-driven", stim_trial_numbers
        )
        noise_sequence_counts = read_data(net, "noise", stim_trial_numbers)
        all_sequence_counts = read_data(net, "all", stim_trial_numbers)

        E_perfect_stimulus_driven_obs = np.mean(
            perfect_stimulus_driven_sequence_counts, axis=(0, 1)
        )
        E_noise_obs = np.mean(noise_sequence_counts, axis=(0, 1))
        E_all_obs = np.mean(all_sequence_counts, axis=(0, 1))
        E_gf_obs = E_all_obs - E_perfect_stimulus_driven_obs - E_noise_obs
        p_gf_obs = 1 - np.exp(-E_gf_obs)

        p = net["connection_prob"]
        p_e = net["ensemble_participation_prob"]
        N = net["ensemble_size"]
        L = net["dendrite_length"]
        M = np.arange(MIN_M, MAX_M + 1)
        delta = net["delta"]
        t_pop_pre = net["population_size"]
        R = net["background_rate"]
        D = net["sequence_time_step"]
        noise_prob = 1 - np.exp(-R * D)
        num_synapses = int(p * t_pop_pre)

        E_perfect_stimulus_driven_seq_exp = (
            p * p_e * N * np.power(p * p_e * N * delta / L, M - 1)
        )
        E_noise_seq_exp = (
            noise_prob
            * (num_synapses - p * N)
            * np.power(noise_prob * (num_synapses - p * N) * delta / L, M - 1)
        )
        E_all_seq_exp = (p * p_e * N + noise_prob * (num_synapses - p * N)) * np.power(
            (p * p_e * N + noise_prob * (num_synapses - p * N)) * delta / L, M - 1
        )
        E_gf_seq_exp = (
            E_all_seq_exp - E_perfect_stimulus_driven_seq_exp - E_noise_seq_exp
        )
        p_gf_seq_exp = 1 - np.exp(-E_gf_seq_exp)
        ax.plot(M, p_gf_seq_exp, ":x", c=colors[i], label=net["name"])
        ax.plot(M, p_gf_obs, "-o", c=colors[i], label=net["name"])

    sns.despine()
    ax.set_yscale("log", nonpositive="mask")
    ax.set_ylim([1e-5, 1e0])
    ax.set_ylabel("p_gap_fill_seq")
    ax.set_xticks([3, 5, 7, 9])
    ax.set_yticklabels([])


def plot_activity_across_trials(net, M, ax, bx):
    """Plot snapshot of true and false positive activity in a set of neurons across trials

    Parameters
    ----------
    ax : axis object
    """
    num_trials = 3
    start_neuron = 20
    start_trial = 4
    num_neurons = 40
    stim_array = np.genfromtxt(
        DATA_PATH + f"stimulus-{net['name']}-run-001.csv", delimiter=","
    )
    connectivity_trial = np.where(stim_array < 0)[0] + 1
    perfect_stimulus_driven_connected_data = read_data(
        net, "perfect-stimulus-driven", connectivity_trial
    )

    stim_trial_numbers = np.where(stim_array > 0)[0] + 1
    all_activity_data = read_data(net, "all", stim_trial_numbers)

    inter_trial_xSpacing = 1
    for t in range(num_trials):
        start_bound = t * (5 + inter_trial_xSpacing)

        # not connected, not active
        X = np.tile(np.arange(0, 5), 8) + start_bound
        Y = np.repeat(np.arange(0, 8), 5)
        X = X[
            np.logical_and(
                perfect_stimulus_driven_connected_data[
                    0, start_neuron : start_neuron + num_neurons, M - MIN_M
                ]
                == 0,
                all_activity_data[
                    t + start_trial,
                    start_neuron : start_neuron + num_neurons,
                    M - MIN_M,
                ]
                == 0,
            )
        ]
        Y = Y[
            np.logical_and(
                perfect_stimulus_driven_connected_data[
                    0, start_neuron : start_neuron + num_neurons, M - MIN_M
                ]
                == 0,
                all_activity_data[
                    t + start_trial,
                    start_neuron : start_neuron + num_neurons,
                    M - MIN_M,
                ]
                == 0,
            )
        ]
        ax.scatter(X, Y, s=30, facecolors="none", edgecolors=[0, 0, 0, 0.5])

        # not connected, active
        X = np.tile(np.arange(0, 5), 8) + start_bound
        Y = np.repeat(np.arange(0, 8), 5)
        X = X[
            np.logical_and(
                perfect_stimulus_driven_connected_data[
                    0, start_neuron : start_neuron + num_neurons, M - MIN_M
                ]
                == 0,
                all_activity_data[
                    t + start_trial,
                    start_neuron : start_neuron + num_neurons,
                    M - MIN_M,
                ]
                > 0,
            )
        ]
        Y = Y[
            np.logical_and(
                perfect_stimulus_driven_connected_data[
                    0, start_neuron : start_neuron + num_neurons, M - MIN_M
                ]
                == 0,
                all_activity_data[
                    t + start_trial,
                    start_neuron : start_neuron + num_neurons,
                    M - MIN_M,
                ]
                > 0,
            )
        ]
        ax.scatter(
            X,
            Y,
            s=30,
            facecolors=[0, 0, 0, 0.5],
            edgecolors=[0, 0, 0, 0.5],
        )

        # connected, not active
        X = np.tile(np.arange(0, 5), 8) + start_bound
        Y = np.repeat(np.arange(0, 8), 5)
        X = X[
            np.logical_and(
                perfect_stimulus_driven_connected_data[
                    0, start_neuron : start_neuron + num_neurons, M - MIN_M
                ]
                > 0,
                all_activity_data[
                    t + start_trial,
                    start_neuron : start_neuron + num_neurons,
                    M - MIN_M,
                ]
                == 0,
            )
        ]
        Y = Y[
            np.logical_and(
                perfect_stimulus_driven_connected_data[
                    0, start_neuron : start_neuron + num_neurons, M - MIN_M
                ]
                > 0,
                all_activity_data[
                    t + start_trial,
                    start_neuron : start_neuron + num_neurons,
                    M - MIN_M,
                ]
                == 0,
            )
        ]
        ax.scatter(X, Y, marker="s", s=30, facecolors="none", edgecolors="r")

        # connected, active
        X = np.tile(np.arange(0, 5), 8) + start_bound
        Y = np.repeat(np.arange(0, 8), 5)
        X = X[
            np.logical_and(
                perfect_stimulus_driven_connected_data[
                    0, start_neuron : start_neuron + num_neurons, M - MIN_M
                ]
                > 0,
                all_activity_data[
                    t + start_trial,
                    start_neuron : start_neuron + num_neurons,
                    M - MIN_M,
                ]
                > 0,
            )
        ]
        Y = Y[
            np.logical_and(
                perfect_stimulus_driven_connected_data[
                    0, start_neuron : start_neuron + num_neurons, M - MIN_M
                ]
                > 0,
                all_activity_data[
                    t + start_trial,
                    start_neuron : start_neuron + num_neurons,
                    M - MIN_M,
                ]
                > 0,
            )
        ]
        ax.scatter(X, Y, marker="s", s=30, facecolors="r", edgecolors="r")
        ax.add_patch(
            mpatches.Rectangle(
                (start_bound - 0.5, -0.5), 5, 8, fc="none", ec=[0, 0, 0, 0.5], lw=1
            )
        )
        ax.text(
            0.16 + t * 0.33,
            -0.05,
            "Trial #%d" % (t + 1),
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        ax.set_xlabel("Trial %d" % (t + 1))
        ax.set_xticks([], [])
        ax.set_yticks([], [])
        ax.spines["right"].set_visible(True)
        ax.spines["top"].set_visible(True)
    ax.set_title(f"{net['name']} M={M}")
    ax.axis("off")

    # plot reliablity of neurons
    reliability = np.mean(all_activity_data > 0, axis=0)
    bins = np.linspace(0, 1, 11)
    perfect_stimulus_driven_connected_neurons = (
        perfect_stimulus_driven_connected_data[0, :, M - MIN_M] > 0
    )
    bx.hist(
        [
            reliability[perfect_stimulus_driven_connected_neurons, M - MIN_M],
            reliability[
                np.logical_not(perfect_stimulus_driven_connected_neurons), M - MIN_M
            ],
        ],
        bins=bins,
        log=True,
        color=[[1, 0, 0, 1], [0, 0, 0, 0.5]],
        label=["connected", "not connected"],
    )
    bx.set_xlabel("Reliablity")
    bx.set_ylabel("#Neurons")
    bx.legend()
    sns.despine()


def plot_activity_stim_vs_background(net, M, ax, bx):
    """Plot snapshot of activity of a set of neurons in stimulus vs background trials

    Parameters
    ----------
    net : dict
        Dictionary with network configuration parameters
    ax : axis object
    """
    num_trials = 2
    start_neuron = 20
    num_neurons = 60
    stim_array = np.genfromtxt(
        DATA_PATH + f"stimulus-{net['name']}-run-001.csv", delimiter=","
    )
    connectivity_trial = np.where(stim_array < 0)[0] + 1
    perfect_stimulus_driven_connected_data = read_data(
        net, "perfect-stimulus-driven", connectivity_trial
    )

    stim_trial_numbers = np.where(stim_array > 0)[0] + 1
    background_trial_numbers = np.where(stim_array == 0)[0] + 1

    inter_trial_xSpacing = 1
    for t in range(num_trials):
        if t == 0:
            # stim trial
            activity_data = read_data(net, "all", [stim_trial_numbers[4]])
            col = "r"
            label = "Stimulus trial"
        else:
            # background trial
            activity_data = read_data(net, "noise", [background_trial_numbers[4]])
            col = [0, 0, 0, 0.5]
            label = "Background trial"

        start_bound = t * (6 + inter_trial_xSpacing)

        # not connected, not active
        X = np.tile(np.arange(0, 6), 10) + start_bound
        Y = np.repeat(np.arange(0, 10), 6)
        X = X[
            np.logical_and(
                perfect_stimulus_driven_connected_data[
                    0, start_neuron : start_neuron + num_neurons, M - MIN_M
                ]
                == 0,
                activity_data[0, start_neuron : start_neuron + num_neurons, M - MIN_M]
                == 0,
            )
        ]
        Y = Y[
            np.logical_and(
                perfect_stimulus_driven_connected_data[
                    0, start_neuron : start_neuron + num_neurons, M - MIN_M
                ]
                == 0,
                activity_data[0, start_neuron : start_neuron + num_neurons, M - MIN_M]
                == 0,
            )
        ]
        ax.scatter(X, Y, s=30, facecolors="none", ec=[0, 0, 0, 0.5])

        # not connected, active
        X = np.tile(np.arange(0, 6), 10) + start_bound
        Y = np.repeat(np.arange(0, 10), 6)
        X = X[
            np.logical_and(
                perfect_stimulus_driven_connected_data[
                    0, start_neuron : start_neuron + num_neurons, M - MIN_M
                ]
                == 0,
                activity_data[0, start_neuron : start_neuron + num_neurons, M - MIN_M]
                > 0,
            )
        ]
        Y = Y[
            np.logical_and(
                perfect_stimulus_driven_connected_data[
                    0, start_neuron : start_neuron + num_neurons, M - MIN_M
                ]
                == 0,
                activity_data[0, start_neuron : start_neuron + num_neurons, M - MIN_M]
                > 0,
            )
        ]
        ax.scatter(X, Y, s=30, fc=[0, 0, 0, 0.5], ec=[0, 0, 0, 0.5])

        # connected, not active
        X = np.tile(np.arange(0, 6), 10) + start_bound
        Y = np.repeat(np.arange(0, 10), 6)
        X = X[
            np.logical_and(
                perfect_stimulus_driven_connected_data[
                    0, start_neuron : start_neuron + num_neurons, M - MIN_M
                ]
                > 0,
                activity_data[0, start_neuron : start_neuron + num_neurons, M - MIN_M]
                == 0,
            )
        ]
        Y = Y[
            np.logical_and(
                perfect_stimulus_driven_connected_data[
                    0, start_neuron : start_neuron + num_neurons, M - MIN_M
                ]
                > 0,
                activity_data[0, start_neuron : start_neuron + num_neurons, M - MIN_M]
                == 0,
            )
        ]
        ax.scatter(X, Y, marker="s", s=30, facecolors="none", edgecolors=col)

        # connected, active
        X = np.tile(np.arange(0, 6), 10) + start_bound
        Y = np.repeat(np.arange(0, 10), 6)
        X = X[
            np.logical_and(
                perfect_stimulus_driven_connected_data[
                    0, start_neuron : start_neuron + num_neurons, M - MIN_M
                ]
                > 0,
                activity_data[0, start_neuron : start_neuron + num_neurons, M - MIN_M]
                > 0,
            )
        ]
        Y = Y[
            np.logical_and(
                perfect_stimulus_driven_connected_data[
                    0, start_neuron : start_neuron + num_neurons, M - MIN_M
                ]
                > 0,
                activity_data[0, start_neuron : start_neuron + num_neurons, M - MIN_M]
                > 0,
            )
        ]
        ax.scatter(X, Y, marker="s", s=30, facecolors=col, edgecolors=col)
        ax.add_patch(
            mpatches.Rectangle(
                (start_bound - 0.5, -0.5), 6, 10, fc="none", ec=[0, 0, 0, 0.5], lw=1
            )
        )
        ax.text(
            0.25 + t * 0.5,
            -0.05,
            f"{label}",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        ax.set_xticks([], [])
        ax.set_yticks([], [])
        ax.spines["right"].set_visible(True)
        ax.spines["top"].set_visible(True)
    ax.set_title(f"{net['name']} M={M}")
    ax.axis("off")


def plot_ratio_perfect_stimulus_driven_active_seq_to_any_seq(ax):
    """Plot probability ratio of perfect stimulus driven active sequences to any sequence

    Parameters
    ----------
    ax : axis object
    """
    for i, net in enumerate(networks):
        stim_array = np.genfromtxt(
            DATA_PATH + f"stimulus-{net['name']}-run-001.csv", delimiter=","
        )
        stim_trial_numbers = np.where(stim_array > 0)[0] + 1
        perfect_stimulus_driven_sequence_counts = read_data(
            net, "perfect-stimulus-driven", stim_trial_numbers
        )
        p_obs_aposs = np.mean(perfect_stimulus_driven_sequence_counts > 0, axis=(0, 1))

        all_sequence_counts = read_data(net, "all", stim_trial_numbers)
        p_obs_any = np.mean(all_sequence_counts > 0, axis=(0, 1))

        p = net["connection_prob"]
        p_e = net["ensemble_participation_prob"]
        N = net["ensemble_size"]
        L = net["dendrite_length"]
        M = np.arange(MIN_M, MAX_M + 1)
        delta = net["delta"]
        t_pop_pre = net["population_size"]
        R = net["background_rate"]
        D = net["sequence_time_step"]
        noise_prob = 1 - np.exp(-R * D)
        num_synapses = int(p * t_pop_pre)

        E_seq_exp_aposs = p * p_e * N * np.power(p * p_e * N * delta / L, M - 1)
        p_seq_exp_aposs = 1 - np.exp(-E_seq_exp_aposs)

        E_seq_exp_any = (p * p_e * N + noise_prob * (num_synapses - p * N)) * np.power(
            (p * p_e * N + noise_prob * (num_synapses - p * N)) * delta / L, M - 1
        )
        p_seq_exp_any = 1 - np.exp(-E_seq_exp_any)
        ax.plot(
            M,
            np.divide(p_seq_exp_aposs, p_seq_exp_any),
            ":x",
            c=colors[i],
            label=net["name"],
        )
        ax.plot(
            M, np.divide(p_obs_aposs, p_obs_any), "-o", c=colors[i], label=net["name"]
        )
    sns.despine()
    ax.set_yscale("log", nonpositive="mask")
    ax.set_ylim([1e-5, 1e0])
    ax.set_ylabel("p_aposs:p_any_seq")
    ax.set_xticks([3, 5, 7, 9])
    ax.set_yticklabels([])
    ax.set_xlabel("Sequence length M")


def plot_legend(ax):
    """Plot legend

    Parameters
    ----------
    ax : axis object
    """
    ax.text(
        0.4,
        0.9,
        "Legend for D,E,F,G,H,I",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    legend_handles = []
    eqn_line = mlines.Line2D(
        [], [], marker="x", ls=":", color="k", label="Approximation"
    )
    legend_handles.append(eqn_line)
    sim_line = mlines.Line2D(
        [],
        [],
        marker="o",
        color="k",
        label="Simulation",
    )
    legend_handles.append(sim_line)
    for i, net in enumerate(networks):
        patch = mpatches.Patch(color=colors[i], label=net["name"])
        legend_handles.append(patch)
    ax.legend(handles=legend_handles, loc="center", frameon=False)
    ax.axis("off")


def main():
    figure = plt.figure(figsize=(10, 10), constrained_layout=True)
    spec = figure.add_gridspec(4, 15)
    axa = figure.add_subplot(spec[0, 0:5])
    axb = figure.add_subplot(spec[0, 5:10])
    axc = figure.add_subplot(spec[0, 10:15])
    axd = figure.add_subplot(spec[1, 0:5])
    axe = figure.add_subplot(spec[1, 5:15])
    # axf = figure.add_subplot(spec[1, 10:15])
    axg = figure.add_subplot(spec[2, 0:5])
    axh = figure.add_subplot(spec[2, 5:10])
    axi = figure.add_subplot(spec[2, 10:15])
    axj = figure.add_subplot(spec[3, 0:5])
    axk = figure.add_subplot(spec[3, 5:10])
    axl = figure.add_subplot(spec[3, 10:15])

    plot_sequence_types(axa)
    plot_active_vs_inactive_sequences(axb)
    plot_legend(axc)
    plot_sequences_venn_diagram(axd, axe)
    plot_perfect_stimulus_driven_conn_sequences(axg)
    plot_noise_sequences(axh)
    plot_gap_fill_sequences(axi)
    plot_perfect_stimulus_driven_active_sequences(axj)
    plot_any_sequences(axk)
    plot_ratio_perfect_stimulus_driven_active_seq_to_any_seq(axl)

    plt.figtext(0.05, 0.98, "A", fontsize=16, weight="bold")
    plt.figtext(0.37, 0.98, "B", fontsize=16, weight="bold")
    plt.figtext(0.03, 0.75, "C", fontsize=16, weight="bold")
    plt.figtext(0.01, 0.52, "D", fontsize=16, weight="bold")
    plt.figtext(0.35, 0.52, "E", fontsize=16, weight="bold")
    plt.figtext(0.69, 0.52, "F", fontsize=16, weight="bold")
    plt.figtext(0.01, 0.25, "G", fontsize=16, weight="bold")
    plt.figtext(0.35, 0.25, "H", fontsize=16, weight="bold")
    plt.figtext(0.69, 0.25, "I", fontsize=16, weight="bold")

    plt.savefig("Figure-4.png", bbox_inches="tight")
    plt.savefig("Figure-4.svg", bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
