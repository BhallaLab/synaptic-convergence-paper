#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : plot_figure3.py
# Author            : Bhanu Priya S <bhanups@ncbs.res.in>
# Date              : 09.11.2023
# Last Modified Date: 11.01.2024
# Last Modified By  : Bhanu Priya S <bhanups@ncbs.res.in>

"""Plot Figure-2
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import seaborn as sns
import pandas as pd
from scipy.stats import poisson
from scipy.linalg import expm
from scipy.optimize import curve_fit
import matplotlib.image as mpimg

sys.path.append("../../codes/groups/")
import group_params as prm

sns.set(style="ticks")
sns.set_context("paper")

DATA_PATH = "../../codes/groups/data/"
CAM_MODEL_DATA_PATH = "../../codes/groups/cam/data/"

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
MIN_POWER = prm.MIN_POWER
MAX_POWER = 3

T_POP_POST = 4000
NUM_RUNS = 100


def plot_group_types(ax):
    """Display different types of groups from ensembles

    Parameters
    ----------
    ax : axis object
    """
    ax.imshow(mpimg.imread("stimulus_driven_groups.png"))
    ax.text(
        -0.05,
        0.13,
        "Partially-mixed\ngroup",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.text(
        -0.06,
        0.74,
        "Fully-mixed\ngroup",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.text(
        0.95,
        0.19,
        "Homogeneous\ngroup",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.set_title("Stimulus-driven groups")
    ax.axis("off")


def plot_active_vs_inactive_groups(ax):
    """Distinguishing active groups from inactive groups

    Parameters
    ----------
    ax : axis object
    """
    ax.imshow(mpimg.imread("active_vs_inactive_groups_and_noise_groups.png"))
    ax.text(
        -0.04,
        0.13,
        "Active\ngroup",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.text(
        0.95,
        0.19,
        "Inactive\ngroup",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.text(
        -0.03,
        0.74,
        "Noise\ngroup",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.axis("off")


def plot_groups_venn_diagram(ax, bx):
    """Distinguishing the venn diagram of groups

    Parameters
    ----------
    ax : axis object
    bx : axis object
    """
    ax.imshow(mpimg.imread("groups_venn_diagram.png"))

    ax.plot([0.17, 0.11], [0.88, 0.91], color="k", transform=ax.transAxes)
    ax.text(
        0.08,
        0.92,
        "U",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.plot([0.71, 0.79], [0.68, 0.73], color="k", transform=ax.transAxes)
    ax.text(
        0.85,
        0.75,
        "CSD",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.plot([0.66, 0.79], [0.40, 0.28], color="k", transform=ax.transAxes)
    ax.text(
        0.85,
        0.25,
        "ASD",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.plot([0.31, 0.13], [0.50, 0.50], color="k", transform=ax.transAxes)
    ax.text(
        0.07,
        0.50,
        "CFM",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.plot([0.41, 0.41], [0.26, 0.48], color="k", transform=ax.transAxes)
    ax.text(
        0.41,
        0.22,
        "AFM",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.axis("off")
    legend_handles = []
    neuron_types_colors = [
        "#000000FF",
        "#44AA99FF",
        "#AA4499FF",
        "#8877DEFF",
        "#88CCEEFF",
    ]
    neuron_types = [
        "$U$ - All postsynaptic neurons",
        "$CSD$ - Neurons receiving grouped\nconnections from stimulus-driven neurons",
        "$ASD$ - Neurons receiving active grouped\ninputs from stimulus-driven neurons",
        "$CFM$ - Neurons receiving\nfully-mixed grouped connectivity",
        "$AFM$ - Neurons receiving active\nfully-mixed grouped inputs",
    ]
    for i, nt in enumerate(neuron_types):
        if i == 4:
            patch = mpatches.Patch(
                edgecolor="none", facecolor=neuron_types_colors[i], label=nt
            )
        else:
            patch = mpatches.Patch(
                edgecolor=neuron_types_colors[i], facecolor="none", label=nt
            )
        legend_handles.append(patch)
    bx.legend(handles=legend_handles, loc="center", frameon=False)
    bx.axis("off")


def plot_ca_cam_model(ax):
    """Display Ca-CaM reaction scheme

    Parameters
    ----------
    ax : axis object
    """
    ax.imshow(mpimg.imread("ca_cam_model.png"))
    ax.set_title("Ca-CaM Model")
    ax.axis("off")
    ax.axis("off")


def plot_ca_cam_traces(ax):
    """Plot Ca-CaM traces

    Parameters
    ----------
    ax : axis object
    """
    grouped_ca_trace = pd.read_csv(
        CAM_MODEL_DATA_PATH + "sample_grouped_output_seq_dx_0.000002_seq_dt_0.000.csv",
        header=None,
        sep=",",
    )

    dispersed_ca_trace = pd.read_csv(
        CAM_MODEL_DATA_PATH
        + "sample_dispersed_output_seq_dx_0.000010_seq_dt_0.000.csv",
        header=None,
        sep=",",
    )

    start_ind = int(grouped_ca_trace.shape[1] * 0.2)
    end_ind = int(grouped_ca_trace.shape[1] * 0.6)

    ax.plot(
        np.arange(start_ind, end_ind) * 0.001,
        grouped_ca_trace.iloc[0][start_ind:end_ind] * 1000,
        label="grouped\n(dx=2\u03bcm)",
    )
    ax.plot(
        np.arange(start_ind, end_ind) * 0.001,
        dispersed_ca_trace.iloc[0][start_ind:end_ind] * 1000,
        label="dispersed\n(dx=10\u03bcm)",
    )
    ax.legend(frameon=False)
    ax.set_xlabel("Time(s)")
    ax.set_ylabel("Ca4_CaM conc (\u03bcM)")


def read_data(net, prefix, trial_numbers):
    """Read group counts data across trials

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
        3D array of group counts of shape num_trials x (t_pop_post x num_runs) x num_group_sizes
    """
    data_matrix = np.zeros(
        (len(trial_numbers), T_POP_POST * NUM_RUNS, MAX_M - MIN_M + 1), dtype=np.int32
    )
    for t, t_num in enumerate(trial_numbers):
        for r in range(NUM_RUNS):
            filename = DATA_PATH + "%s-groups-%s-run-%03d-trial-%03d.csv" % (
                prefix,
                net["name"],
                r + 1,
                t_num,
            )
            data = np.asarray(pd.read_csv(filename, header=None))
            data_matrix[t, r * T_POP_POST : (r + 1) * T_POP_POST, :] = data

    return data_matrix


def read_stimulus_driven_group_splits_data(net, M):
    """Read stimulus-driven group splits data for group size M

    Parameters
    ----------
    net : dict
        Dictionary of network parameters
    M : int
        Group size

    Returns
    -------
    data_matrix : np.ndarray
        2D array of stimulus-driven group split counts of shape (t_pop_post x num_runs) x num_ensembles
    """
    data_matrix = np.zeros(
        (T_POP_POST * NUM_RUNS, M),
        dtype=np.int32,
    )
    for r in range(NUM_RUNS):
        filename = DATA_PATH + "stimulus-driven-group-splits-%s-run-%03d-M-%d.csv" % (
            net["name"],
            r + 1,
            M,
        )
        data = np.asarray(pd.read_csv(filename, header=None))
        data_matrix[r * T_POP_POST : (r + 1) * T_POP_POST, :] = data

    return data_matrix


def plot_fully_mixed_conn_groups(ax):
    """Plot probability of fully-mixed groups based on connectivity

    Parameters
    ----------
    ax : axis object
    """
    for i, net in enumerate(networks):
        stim_array = np.genfromtxt(
            DATA_PATH + f"stimulus-{net['name']}-run-001.csv", delimiter=","
        )
        stim_trial_numbers = np.where(stim_array < 0)[0] + 1
        fully_mixed_group_counts = read_data(net, "fully-mixed", stim_trial_numbers)
        p_obs = np.mean(fully_mixed_group_counts > 0, axis=(0, 1))

        p = net["connection_prob"]
        N = net["ensemble_size"]
        L = net["dendrite_length"]
        M = np.arange(MIN_M, MAX_M + 1)
        Z = net["zone_length"]
        t_pop_pre = net["population_size"]
        sigma = L / int(p * t_pop_pre)

        def func(x, c):
            out = []
            for m in x:
                out.append(
                    np.log(
                        1
                        - np.power(
                            (1 - (np.power(1 - np.exp(-p * N * Z / L), m))),
                            c * int((L / Z)),
                        )
                    )
                )
            return out

        popt, pcov = curve_fit(func, M[p_obs > 1e-5], np.log(p_obs[p_obs > 1e-5]))

        egrp_exp_Ei = p * N * Z / L
        p_grp_exp = np.zeros(len(M))
        p_grp_exp = np.power(1 - np.exp(-egrp_exp_Ei), M)
        p_total_exp_1 = 1 - np.power((1 - p_grp_exp), int(L / Z))
        p_total_exp_2 = 1 - np.power((1 - p_grp_exp), int(L / sigma))
        ax.fill_between(
            M[p_total_exp_1 > 0],
            p_total_exp_1[p_total_exp_1 > 0],
            p_total_exp_2[p_total_exp_1 > 0],
            color=colors[i],
            alpha=0.3,
            label=net["name"],
        )
        ax.plot(M, p_obs, "-o", c=colors[i], label=net["name"])
        ax.plot(M, np.exp(func(M, *popt)), ":x", c=colors[i], label=net["name"])
        print("cFMG", net["name"], *popt)

    sns.despine()
    ax.set_yscale("log", nonpositive="mask")
    ax.set_ylim([1e-5, 1e0])
    ax.set_ylabel("p_fully_mixed_group")
    ax.set_xticks([3, 5, 7, 9])
    ax.set_xlabel("Group size >=M")


def plot_stimulus_driven_conn_groups(ax):
    """Plot probability of stimulus-driven groups based on connectivity

    Parameters
    ----------
    ax : axis object
    """
    for i, net in enumerate(networks):
        stim_array = np.genfromtxt(
            DATA_PATH + f"stimulus-{net['name']}-run-001.csv", delimiter=","
        )
        stim_trial_numbers = np.where(stim_array < 0)[0] + 1
        stimulus_driven_group_counts = read_data(
            net, "stimulus-driven", stim_trial_numbers
        )
        p_obs = np.mean(stimulus_driven_group_counts > 0, axis=(0, 1))

        p = net["connection_prob"]
        N = net["ensemble_size"]
        L = net["dendrite_length"]
        M = np.arange(MIN_M, MAX_M + 1)
        Z = net["zone_length"]
        t_pop_pre = net["population_size"]
        sigma = L / int(p * t_pop_pre)

        def func(x, c):
            out = []
            for m in x:
                out.append(
                    np.log(
                        1
                        - np.power(
                            1
                            - (
                                1
                                - np.sum(
                                    [
                                        poisson.pmf(k=j, mu=(p * N * m * Z / L))
                                        for j in np.arange(0, m)
                                    ]
                                )
                            ),
                            c * int((L / Z)),
                        )
                    )
                )
            return out

        popt, pcov = curve_fit(func, M[p_obs > 5e-5], np.log(p_obs[p_obs > 5e-5]))

        p_alt_exp = np.zeros(len(M))
        for m in M:
            ealt_all_exp_EM = p * N * m * Z / L
            p_alt_exp[m - MIN_M] = 1 - np.sum(
                [poisson.pmf(k=j, mu=ealt_all_exp_EM) for j in range(0, m)]
            )
        p_total_exp_1 = 1 - np.power((1 - p_alt_exp), int(L / Z))
        p_total_exp_2 = 1 - np.power((1 - p_alt_exp), int(L / sigma))
        ax.fill_between(
            M[p_total_exp_1 > 0],
            p_total_exp_1[p_total_exp_1 > 0],
            p_total_exp_2[p_total_exp_1 > 0],
            color=colors[i],
            alpha=0.3,
            label=net["name"],
        )
        ax.plot(M, p_obs, "-o", c=colors[i], label=net["name"])
        ax.plot(M, np.exp(func(M, *popt)), ":x", c=colors[i], label=net["name"])
        print("cSDG", net["name"], *popt)

    sns.despine()
    ax.set_yscale("log", nonpositive="mask")
    ax.set_ylim([1e-5, 1e0])
    ax.set_ylabel("p_stimulus_driven_group")
    ax.set_xticks([3, 5, 7, 9])
    ax.set_yticklabels([])
    ax.set_xlabel("Group size >=M")


def plot_fully_mixed_active_groups(ax):
    """Plot probability of fully-mixed groups based on activity

    Parameters
    ----------
    ax : axis object
    """
    for i, net in enumerate(networks):
        stim_array = np.genfromtxt(
            DATA_PATH + f"stimulus-{net['name']}-run-001.csv", delimiter=","
        )
        stim_trial_numbers = np.where(stim_array > 0)[0] + 1
        fully_mixed_group_counts = read_data(net, "fully-mixed", stim_trial_numbers)
        p_obs = np.mean(fully_mixed_group_counts > 0, axis=(0, 1))

        p = net["connection_prob"]
        p_e = net["ensemble_participation_prob"]
        N = net["ensemble_size"]
        L = net["dendrite_length"]
        M = np.arange(MIN_M, MAX_M + 1)
        Z = net["zone_length"]
        t_pop_pre = net["population_size"]
        sigma = L / int(p * t_pop_pre)

        def func(x, c):
            out = []
            for m in x:
                out.append(
                    np.log(
                        1
                        - np.power(
                            (1 - (np.power(1 - np.exp(-p * p_e * N * Z / L), m))),
                            c * int((L / Z)),
                        )
                    )
                )
            return out

        popt, pcov = curve_fit(func, M[p_obs > 1e-5], np.log(p_obs[p_obs > 1e-5]))

        egrp_exp_Ei = p * p_e * N * Z / L
        p_grp_exp = np.zeros(len(M))
        p_grp_exp = np.power(1 - np.exp(-egrp_exp_Ei), M)
        p_total_exp_1 = 1 - np.power((1 - p_grp_exp), int(L / Z))
        p_total_exp_2 = 1 - np.power((1 - p_grp_exp), int(L / sigma))
        ax.fill_between(
            M[p_total_exp_1 > 0],
            p_total_exp_1[p_total_exp_1 > 0],
            p_total_exp_2[p_total_exp_1 > 0],
            color=colors[i],
            alpha=0.3,
            label=net["name"],
        )
        p_total_exp_fit = []
        ax.plot(M, p_obs, "-o", c=colors[i], label=net["name"])
        ax.plot(M, np.exp(func(M, *popt)), ":x", c=colors[i], label=net["name"])
        print("aFMG", net["name"], *popt)

    sns.despine()
    ax.set_yscale("log", nonpositive="mask")
    ax.set_ylim([1e-5, 1e0])
    ax.set_ylabel("p_fully_mixed_group")
    ax.set_xticks([3, 5, 7, 9])
    ax.set_xlabel("Group size >=M")


def plot_stimulus_driven_active_groups(ax):
    """Plot probability of stimulus-driven groups based on activity

    Parameters
    ----------
    ax : axis object
    """
    for i, net in enumerate(networks):
        stim_array = np.genfromtxt(
            DATA_PATH + f"stimulus-{net['name']}-run-001.csv", delimiter=","
        )
        stim_trial_numbers = np.where(stim_array > 0)[0] + 1
        stimulus_driven_group_counts = read_data(
            net, "stimulus-driven", stim_trial_numbers
        )
        p_obs = np.mean(stimulus_driven_group_counts > 0, axis=(0, 1))

        p = net["connection_prob"]
        p_e = net["ensemble_participation_prob"]
        N = net["ensemble_size"]
        L = net["dendrite_length"]
        M = np.arange(MIN_M, MAX_M + 1)
        Z = net["zone_length"]
        t_pop_pre = net["population_size"]
        sigma = L / int(p * t_pop_pre)

        def func(x, c):
            out = []
            for m in x:
                out.append(
                    np.log(
                        1
                        - np.power(
                            1
                            - (
                                1
                                - np.sum(
                                    [
                                        poisson.pmf(k=j, mu=(p * p_e * N * m * Z / L))
                                        for j in np.arange(0, m)
                                    ]
                                )
                            ),
                            c * int((L / Z)),
                        )
                    )
                )
            return out

        popt, pcov = curve_fit(func, M[p_obs > 5e-5], np.log(p_obs[p_obs > 5e-5]))

        p_alt_exp = np.zeros(len(M))
        for m in M:
            ealt_all_exp_EM = p * p_e * N * m * Z / L
            p_alt_exp[m - MIN_M] = 1 - np.sum(
                [poisson.pmf(k=i, mu=ealt_all_exp_EM) for i in range(0, m)]
            )
        p_total_exp_1 = 1 - np.power((1 - p_alt_exp), int(L / Z))
        p_total_exp_2 = 1 - np.power((1 - p_alt_exp), int(L / sigma))
        ax.fill_between(
            M[p_total_exp_1 > 0],
            p_total_exp_1[p_total_exp_1 > 0],
            p_total_exp_2[p_total_exp_1 > 0],
            color=colors[i],
            alpha=0.3,
            label=net["name"],
        )
        p_total_exp_fit = []
        ax.plot(M, p_obs, "-o", c=colors[i], label=net["name"])
        ax.plot(M, np.exp(func(M, *popt)), ":x", c=colors[i], label=net["name"])
        print("aSDG", net["name"], *popt)

    sns.despine()
    ax.set_yscale("log", nonpositive="mask")
    ax.set_ylim([1e-5, 1e0])
    ax.set_ylabel("p_stimulus_driven_group")
    ax.set_xticks([3, 5, 7, 9])
    ax.set_yticklabels([])
    ax.set_xlabel("Group size >=M")


def plot_noise_groups(ax):
    """Plot probability of noise groups based on activity

    Parameters
    ----------
    ax : axis object
    """
    for i, net in enumerate(networks):
        stim_array = np.genfromtxt(
            DATA_PATH + f"stimulus-{net['name']}-run-001.csv", delimiter=","
        )
        stim_trial_numbers = np.where(stim_array > 0)[0] + 1
        noise_group_counts = read_data(net, "noise", stim_trial_numbers)
        p_obs = np.mean(noise_group_counts > 0, axis=(0, 1))

        p = net["connection_prob"]
        N = net["ensemble_size"]
        L = net["dendrite_length"]
        M = np.arange(MIN_M, MAX_M + 1)
        Z = net["zone_length"]
        t_pop_pre = net["population_size"]
        sigma = L / int(p * t_pop_pre)
        R = net["background_rate"]
        D = net["input_duration"]

        p_noise_exp = np.zeros(len(M))
        p_bg = 1 - np.exp(-R * D)
        for m in M:
            ENoise_exp = p_bg * Z * ((1 / sigma) - (p * N * m / L))
            p_noise_exp[m - MIN_M] = 1 - np.sum(
                [poisson.pmf(k=i, mu=ENoise_exp) for i in range(0, m)]
            )
        p_total_exp_1 = 1 - np.power((1 - p_noise_exp), int(L / Z))
        p_total_exp_2 = 1 - np.power((1 - p_noise_exp), int(L / sigma))
        ax.fill_between(
            M[p_total_exp_1 > 0],
            p_total_exp_1[p_total_exp_1 > 0],
            p_total_exp_2[p_total_exp_1 > 0],
            color=colors[i],
            alpha=0.3,
            label=net["name"],
        )

        def func(x, c):
            out = []
            for m in x:
                out.append(
                    np.log(
                        1
                        - np.power(
                            1
                            - (
                                1
                                - np.sum(
                                    [
                                        poisson.pmf(
                                            k=j,
                                            mu=(
                                                p_bg
                                                * Z
                                                * ((1 / sigma) - (p * N * m / L))
                                            ),
                                        )
                                        for j in np.arange(0, m)
                                    ]
                                )
                            ),
                            c * int((L / Z)),
                        )
                    )
                )
            return out

        popt, pcov = curve_fit(func, M[p_obs > 5e-5], np.log(p_obs[p_obs > 5e-5]))

        p_total_exp_fit = []
        ax.plot(M, p_obs, "-o", c=colors[i], label=net["name"])
        ax.plot(M, np.exp(func(M, *popt)), ":x", c=colors[i], label=net["name"])
        print("NG", net["name"], *popt)

    sns.despine()
    ax.set_yscale("log", nonpositive="mask")
    ax.set_ylim([1e-5, 1e0])
    ax.set_xlabel("Group size >=M")
    ax.set_ylabel("p_noise_group")
    ax.set_xticks([3, 5, 7, 9])
    ax.set_yticklabels([])


def plot_any_groups(ax):
    """Plot probability of any group based on activity

    Parameters
    ----------
    ax : axis object
    """
    for i, net in enumerate(networks):
        stim_array = np.genfromtxt(
            DATA_PATH + f"stimulus-{net['name']}-run-001.csv", delimiter=","
        )
        stim_trial_numbers = np.where(stim_array > 0)[0] + 1
        all_group_counts = read_data(net, "all", stim_trial_numbers)
        p_obs = np.mean(all_group_counts > 0, axis=(0, 1))

        p = net["connection_prob"]
        p_e = net["ensemble_participation_prob"]
        N = net["ensemble_size"]
        L = net["dendrite_length"]
        M = np.arange(MIN_M, MAX_M + 1)
        Z = net["zone_length"]
        t_pop_pre = net["population_size"]
        sigma = L / int(p * t_pop_pre)
        R = net["background_rate"]
        D = net["input_duration"]

        p_any_exp = np.zeros(len(M))
        p_bg = 1 - np.exp(-R * D)
        for m in M:
            EAll_exp = p_bg * Z * ((1 / sigma) - (p * N * m / L)) + (
                p * p_e * N * m * Z / L
            )
            p_any_exp[m - MIN_M] = 1 - np.sum(
                [poisson.pmf(k=i, mu=EAll_exp) for i in range(0, m)]
            )
        p_any_exp_1 = 1 - np.power((1 - p_any_exp), int(L / Z))
        p_any_exp_2 = 1 - np.power((1 - p_any_exp), int(L / sigma))
        ax.fill_between(
            M[p_any_exp_1 > 0],
            p_any_exp_1[p_any_exp_1 > 0],
            p_any_exp_2[p_any_exp_1 > 0],
            color=colors[i],
            alpha=0.3,
            label=net["name"],
        )

        def func(x, c):
            out = []
            for m in x:
                out.append(
                    np.log(
                        1
                        - np.power(
                            1
                            - (
                                1
                                - np.sum(
                                    [
                                        poisson.pmf(
                                            k=j,
                                            mu=(
                                                p_bg
                                                * Z
                                                * ((1 / sigma) - (p * N * m / L))
                                                + (p * p_e * N * m * Z / L)
                                            ),
                                        )
                                        for j in np.arange(0, m)
                                    ]
                                )
                            ),
                            c * int((L / Z)),
                        )
                    )
                )
            return out

        popt, pcov = curve_fit(func, M[p_obs > 5e-5], np.log(p_obs[p_obs > 5e-5]))

        p_total_exp_fit = []
        ax.plot(M, p_obs, "-o", c=colors[i], label=net["name"])
        ax.plot(M, np.exp(func(M, *popt)), ":x", c=colors[i], label=net["name"])
        print("AG", net["name"], *popt)

    sns.despine()
    ax.set_yscale("log", nonpositive="mask")
    ax.set_ylim([1e-5, 1e0])
    ax.set_xlabel("Group size >=M")
    ax.set_ylabel("p_any_group")
    ax.set_xticks([3, 5, 7, 9])
    ax.set_yticklabels([])


def plot_ratio_stimulus_driven_to_any_group(ax):
    """Plot ratio of probability of stimulus-driven active groups to any groups
    Parameters
    ----------
    ax : axis object
    """

    for i, net in enumerate(networks):
        stim_array = np.genfromtxt(
            DATA_PATH + f"stimulus-{net['name']}-run-001.csv", delimiter=","
        )
        stim_trial_numbers = np.where(stim_array > 0)[0] + 1
        stimulus_driven_group_counts = read_data(
            net, "stimulus-driven", stim_trial_numbers
        )
        p_stimulus_driven_obs = np.mean(stimulus_driven_group_counts > 0, axis=(0, 1))
        all_group_counts = read_data(net, "all", stim_trial_numbers)
        p_any_obs = np.mean(all_group_counts > 0, axis=(0, 1))

        p = net["connection_prob"]
        p_e = net["ensemble_participation_prob"]
        N = net["ensemble_size"]
        L = net["dendrite_length"]
        M = np.arange(MIN_M, MAX_M + 1)
        Z = net["zone_length"]
        t_pop_pre = net["population_size"]
        sigma = L / int(p * t_pop_pre)
        R = net["background_rate"]
        D = net["input_duration"]
        p_bg = 1 - np.exp(-R * D)

        def func_aSDG(x, c):
            out = []
            for m in x:
                out.append(
                    np.log(
                        1
                        - np.power(
                            1
                            - (
                                1
                                - np.sum(
                                    [
                                        poisson.pmf(k=j, mu=(p * p_e * N * m * Z / L))
                                        for j in np.arange(0, m)
                                    ]
                                )
                            ),
                            c * int((L / Z)),
                        )
                    )
                )
            return out

        popt_aSDG, pcov_aSDG = curve_fit(
            func_aSDG,
            M[p_stimulus_driven_obs > 5e-5],
            np.log(p_stimulus_driven_obs[p_stimulus_driven_obs > 5e-5]),
        )
        p_total_exp_fit_stimulus_driven_group = np.exp(func_aSDG(M, *popt_aSDG))

        def func_any(x, c):
            out = []
            for m in x:
                out.append(
                    np.log(
                        1
                        - np.power(
                            1
                            - (
                                1
                                - np.sum(
                                    [
                                        poisson.pmf(
                                            k=j,
                                            mu=(
                                                p_bg
                                                * Z
                                                * ((1 / sigma) - (p * N * m / L))
                                                + (p * p_e * N * m * Z / L)
                                            ),
                                        )
                                        for j in np.arange(0, m)
                                    ]
                                )
                            ),
                            c * int((L / Z)),
                        )
                    )
                )
            return out

        popt_any, pcov_any = curve_fit(
            func_any, M[p_any_obs > 5e-5], np.log(p_any_obs[p_any_obs > 5e-5])
        )
        p_total_exp_fit_any_group = np.exp(func_any(M, *popt_any))

        ax.plot(
            M,
            np.divide(p_total_exp_fit_stimulus_driven_group, p_total_exp_fit_any_group),
            ":x",
            c=colors[i],
            label=net["name"],
        )
        ax.plot(
            M,
            np.divide(p_stimulus_driven_obs, p_any_obs),
            "-o",
            c=colors[i],
            label=net["name"],
        )

    sns.despine()
    ax.set_yscale("log", nonpositive="mask")
    ax.set_ylim([1e-5, 1e0])
    ax.set_xlabel("Group size >=M")
    ax.set_ylabel("p_ASD_group:p_any_group")
    ax.set_xticks([3, 5, 7, 9])


def plot_stimulus_driven_group_splits(M, ax):
    """Plot splits of stimulus-driven groups for network for a group size M

    Parameters
    ----------
    net : dict
        Dictionary with network configuration parameters
    M : int
        Group size
    ax : axis object
    """
    for i, net in enumerate(networks):
        stim_array = np.genfromtxt(
            DATA_PATH + f"stimulus-{net['name']}-run-001.csv", delimiter=","
        )
        stim_trial_numbers = np.where(stim_array < 0)[0] + 1
        stimulus_driven_group_counts = read_data(
            net, "stimulus-driven", stim_trial_numbers
        )
        stimulus_driven_group_counts_m = stimulus_driven_group_counts[0, :, M - MIN_M]
        stimulus_driven_group_splits = read_stimulus_driven_group_splits_data(net, M)
        split_counts = stimulus_driven_group_splits[
            stimulus_driven_group_counts_m > 0, :
        ]
        ax.plot(range(1, M + 1), np.sum(split_counts, axis=0), "-o", color=colors[i])
    ax.set_xlabel("#Unique ensemble inputs")
    ax.set_ylabel("#Groups")
    ax.set_yscale("log", nonpositive="mask")
    sns.despine()

    return


def plot_legend(ax):
    """Plot legend

    Parameters
    ----------
    ax : axis object
    """
    ax.text(
        0.4,
        0.9,
        "Legend for\nF,G,H,I,J,K,L,M",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    legend_handles = []
    eqn_line = mlines.Line2D([], [], marker="x", ls=":", color="k", label="Analytical")
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
    figure = plt.figure(figsize=(12, 12), constrained_layout=True)
    spec = figure.add_gridspec(4, 12)
    axa = figure.add_subplot(spec[0:2, 0:3])
    axb = figure.add_subplot(spec[0, 3:6])
    axc = figure.add_subplot(spec[0, 6:9])
    axd = figure.add_subplot(spec[0, 9:12])
    # axe = figure.add_subplot(spec[1, 0:3])
    axf = figure.add_subplot(spec[1, 3:6])
    axg = figure.add_subplot(spec[1, 6:9])
    axh = figure.add_subplot(spec[1, 9:12])
    axi = figure.add_subplot(spec[2, 0:3])
    axj = figure.add_subplot(spec[2, 3:6])
    axk = figure.add_subplot(spec[2, 6:9])
    axl = figure.add_subplot(spec[2, 9:12])
    axm = figure.add_subplot(spec[3, 0:3])
    axn = figure.add_subplot(spec[3, 3:6])
    axo = figure.add_subplot(spec[3, 6:9])
    axp = figure.add_subplot(spec[3, 9:12])

    plot_ca_cam_model(axa)
    plot_ca_cam_traces(axb)
    plot_group_types(axc)
    plot_active_vs_inactive_groups(axd)
    plot_groups_venn_diagram(axf, axg)
    plot_legend(axh)
    plot_fully_mixed_conn_groups(axi)
    plot_stimulus_driven_conn_groups(axj)
    plot_noise_groups(axk)
    plot_fully_mixed_active_groups(axm)
    plot_stimulus_driven_active_groups(axn)
    plot_any_groups(axl)
    plot_stimulus_driven_group_splits(4, axo)
    plot_ratio_stimulus_driven_to_any_group(axp)

    plt.figtext(0.01, 0.98, "A", fontsize=12, weight="bold")
    plt.figtext(0.25, 0.98, "B", fontsize=12, weight="bold")
    plt.figtext(0.51, 0.98, "C", fontsize=12, weight="bold")
    plt.figtext(0.78, 0.98, "D", fontsize=12, weight="bold")
    plt.figtext(0.27, 0.70, "E", fontsize=12, weight="bold")
    plt.figtext(0.01, 0.50, "F", fontsize=12, weight="bold")
    plt.figtext(0.27, 0.50, "G", fontsize=12, weight="bold")
    plt.figtext(0.50, 0.50, "H", fontsize=12, weight="bold")
    plt.figtext(0.78, 0.50, "I", fontsize=12, weight="bold")
    plt.figtext(0.01, 0.24, "J", fontsize=12, weight="bold")
    plt.figtext(0.27, 0.24, "K", fontsize=12, weight="bold")
    plt.figtext(0.50, 0.24, "L", fontsize=12, weight="bold")
    plt.figtext(0.76, 0.24, "M", fontsize=12, weight="bold")

    # Annotations
    plt.figtext(
        0.25, 0.52, "Connectivity based", fontsize="large", horizontalalignment="center"
    )
    plt.figtext(
        0.25, 0.26, "Activity based", fontsize="large", horizontalalignment="center"
    )
    plt.savefig("Figure-2.png", bbox_inches="tight")
    plt.savefig("Figure-2.svg", bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
