#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : plot_figure3a.py
# Author            : Bhanu Priya S <bhanups@ncbs.res.in>
# Date              : 27.02.2024
# Last Modified Date: 27.02.2024
# Last Modified By  : Bhanu Priya S <bhanups@ncbs.res.in>

"""Plot Figure-3a
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib import cm
from matplotlib.colors import LogNorm
from pyparsing import col
import seaborn as sns
import pandas as pd


sys.path.append("../../codes/groups/")
import group_params as prm

sys.path.append("../../../utils/")
from scalebars import add_scalebar

sns.set(style="ticks")
sns.set_context("paper")

DATA_PATH = "../../codes/groups/activation_data/"

networks = [
    prm.hippo_chem,
    prm.hippo_cicr,
    prm.hippo_elec,
    prm.cortex_chem,
    prm.cortex_cicr,
    prm.cortex_elec,
]

MIN_M = prm.MIN_M
MAX_M = 4

T_POP_POST = 4000
NUM_RUNS = 100


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


def read_activation_data(net, trial_numbers, power, nonlinearity):
    """Read group counts data across trials

    Parameters
    ----------
    net : dict
        Dictionary of network parameters
    trial_numbers : np.array
        Array of trial numbers to read
    power : int
        power of the nonlinearity to be read

    Returns
    -------
    data_matrix : np.ndarray
        3D array of group counts of shape num_trials x (t_pop_post x num_runs) x num_group_sizes
    """
    data_matrix = np.zeros(
        (len(trial_numbers), T_POP_POST * NUM_RUNS, MAX_M - MIN_M + 1),
        dtype=np.double,
    )
    for t, t_num in enumerate(trial_numbers):
        for r in range(NUM_RUNS):
            filename = (
                DATA_PATH
                + "activation-of-groups-%s-run-%03d-trial-%03d-%s-%d.csv"
                % (
                    net["name"],
                    r + 1,
                    t_num,
                    nonlinearity,
                    power,
                )
            )
            data = np.asarray(pd.read_csv(filename, header=None))
            data_matrix[t, r * T_POP_POST : (r + 1) * T_POP_POST, :] = data

    return data_matrix


def weibull_stretched_exponential(x, a, b, p):
    return 1 - np.power(a, -np.power((x / b), p))


def plot_weibull(ax, M, a, powers):
    x = np.linspace(0, M + 2, 100)
    colors = sns.color_palette("muted", n_colors=len(powers))
    for p, power in enumerate(powers):
        y = weibull_stretched_exponential(x, a, M, power)
        ax.plot(x, y, ls="--", color=colors[p])
        ax.plot(
            np.arange(M + 2),
            weibull_stretched_exponential(np.arange(M + 2), a, M, power),
            ls="",
            marker="o",
            color=colors[p],
            label=power,
        )

    sns.despine()
    ax.set_ylabel("Activation")
    ax.set_xticks(range(M + 2))
    ax.set_xlabel("#Inputs within zone length Z")
    ax.set_title("M = 4")
    ax.legend(title="Power p", loc="lower right", frameon=False)
    ax.text(
        0.05,
        0.85,
        r"$1 - {10}^{{-(\frac{x}{M})}^{\rho}}$",
        horizontalalignment="left",
        verticalalignment="center",
        transform=ax.transAxes,
        fontsize="x-large",
    )


def plot_activation_data(net, M, nonlinearity, ax, powers):
    """Plot neuron activation data

    Parameters
    ----------
    net : dict
        Dictionary with network configuration parameters
    ax : axis object
    """
    stim_array = np.genfromtxt(
        DATA_PATH + f"stimulus-{net['name']}-run-001.csv", delimiter=","
    )
    connectivity_trial = np.where(stim_array < 0)[0] + 1
    fully_mixed_connected_data = read_data(net, "fully-mixed", connectivity_trial)[
        0, :, M - MIN_M
    ]
    stimulus_driven_connected_data = read_data(
        net, "stimulus-driven", connectivity_trial
    )[0, :, M - MIN_M]

    stim_trial_numbers = np.where(stim_array > 0)[0] + 1
    any_active_data = read_data(net, "all", stim_trial_numbers)[:, :, M - MIN_M]

    for p, power in enumerate(powers):
        activation_data = read_activation_data(
            net, stim_trial_numbers, power, nonlinearity
        )[:, :, M - MIN_M]

        fully_mixed_group_activation_data = activation_data[
            :, fully_mixed_connected_data > 0
        ]
        stimulus_driven_group_activation_data = activation_data[
            :, stimulus_driven_connected_data > 0
        ]
        non_fully_mixed_group_activation_data = activation_data[
            :, fully_mixed_connected_data == 0
        ]
        non_stimulus_driven_group_activation_data = activation_data[
            :, stimulus_driven_connected_data == 0
        ]

        half_max_actn = np.max(activation_data) * 0.5

        for n in np.arange(non_stimulus_driven_group_activation_data.shape[1])[:3]:
            ax[p].plot(
                np.arange(non_stimulus_driven_group_activation_data.shape[0]),
                n * half_max_actn + non_stimulus_driven_group_activation_data[:, n],
                color="m",
                label="not stimulus-driven connected",
            )
            hit_trials = np.where(
                any_active_data[:, stimulus_driven_connected_data == 0][:, n] > 0
            )[0]
            ax[p].scatter(
                hit_trials,
                n * half_max_actn
                + non_stimulus_driven_group_activation_data[hit_trials, n],
                c="k",
                alpha=0.5,
                s=10,
            )

        for n in np.arange(stimulus_driven_group_activation_data.shape[1])[:3]:
            ax[p].plot(
                np.arange(stimulus_driven_group_activation_data.shape[0]),
                (4 + p + n) * half_max_actn
                + stimulus_driven_group_activation_data[:, n],
                color="c",
                label="stimulus-driven connected",
            )
            hit_trials = np.where(
                any_active_data[:, stimulus_driven_connected_data > 0][:, n] > 0
            )[0]
            ax[p].scatter(
                hit_trials,
                (4 + p + n) * half_max_actn
                + stimulus_driven_group_activation_data[hit_trials, n],
                c="k",
                alpha=0.5,
                s=10,
            )

        ax[p].set_axis_off()
        x_scalebar_size = np.diff(ax[p].xaxis.get_majorticklocs())[0]
        y_scalebar_size = np.diff(ax[p].yaxis.get_majorticklocs())[0]
        add_scalebar(
            ax[p],
            matchx=False,
            matchy=False,
            sizex=x_scalebar_size,
            sizey=y_scalebar_size,
            labelx=f"{x_scalebar_size} Trials",
            labely=f"{y_scalebar_size} a.u.",
            loc="center left",
        )
        ax[p].text(
            -0.04,
            0.5,
            "Activation",
            horizontalalignment="center",
            verticalalignment="center",
            rotation=90,
            transform=ax[p].transAxes,
        )

        bins = np.linspace(
            np.min(np.mean(activation_data, axis=0)),
            np.max(np.mean(activation_data, axis=0)),
            11,
        )

        ax[2 * p + 2].hist(
            [
                np.mean(stimulus_driven_group_activation_data, axis=0),
                np.mean(non_stimulus_driven_group_activation_data, axis=0),
            ],
            bins=bins,
            log=True,
            color=["c", "m"],
            label=["stimulus-driven connected", "not stimulus-driven connected"],
        )
        ax[2 * p + 2].set_xlabel("Trial averaged activation")
        ax[2 * p + 2].axvline(bins[2], ls="--", color="k")
        if p == 0:
            ax[2 * p + 2].set_ylabel("#Neurons")
        else:
            ax[2 * p + 2].set_yticklabels([])

        ax[2 * p + 3].hist(
            [
                np.mean(fully_mixed_group_activation_data, axis=0),
                np.mean(non_fully_mixed_group_activation_data, axis=0),
            ],
            bins=bins,
            log=True,
            color=[[1, 0, 0, 1], [0, 0, 0, 0.5]],
            label=["fully-mixed connected", "not fully-mixed connected"],
        )
        ax[2 * p + 3].set_xlabel("Trial averaged activation")
        ax[2 * p + 3].axvline(bins[2], ls="--", color="k")
        ax[2 * p + 3].set_yticklabels([])

        bins = np.linspace(
            np.min(activation_data),
            np.max(activation_data),
            11,
        )

        for d, data in enumerate(
            [
                stimulus_driven_group_activation_data,
                non_stimulus_driven_group_activation_data,
            ]
        ):
            trial_counts = []
            for n in np.arange(data.shape[1]):
                counts, _ = np.histogram(data[:, n], bins=bins)
                trial_counts.append(counts)
            trial_counts = np.array(trial_counts)
            neuron_count_list = []
            for i in range(trial_counts.shape[1]):
                cell_counts, _ = np.histogram(
                    trial_counts[:, i],
                    bins=np.arange(0.5, activation_data.shape[0] + 0.5),
                )
                neuron_count_list.append(cell_counts)
            neuron_count_matrix = np.array(neuron_count_list).T

            yticks = np.arange(9.5, activation_data.shape[0], 10)
            yticklabels = [
                f"{i:0.3f}" for i in np.arange(10 / activation_data.shape[0], 1, 0.1)
            ]
            xticks = np.arange(0, len(bins), 2)
            xticklabels = [f"{i:.2f}" for i in bins[::2]]

            if p == 0 and d == 0:
                heatmap = sns.heatmap(
                    neuron_count_matrix,
                    ax=ax[2 * p + 6 + d],
                    norm=LogNorm(1e0, activation_data.shape[1]),
                    cmap=cm.viridis,
                    cbar_kws={
                        "label": "#Neurons",
                        "ticks": [1e5, 1e3, 1e1],
                    },
                    cbar_ax=ax[2 * p + 10],
                    yticklabels=False,
                    xticklabels=False,
                )
            else:
                sns.heatmap(
                    neuron_count_matrix,
                    ax=ax[2 * p + 6 + d],
                    norm=LogNorm(1e0, activation_data.shape[1]),
                    cmap=cm.viridis,
                    cbar=False,
                    yticklabels=False,
                    xticklabels=False,
                )
            ax[2 * p + 6 + d].invert_yaxis()
            ax[2 * p + 6 + d].set_xlabel("Activation")
            ax[2 * p + 6 + d].set_yticks(yticks)
            ax[2 * p + 6 + d].set_xticks(xticks)
            ax[2 * p + 6 + d].set_xticklabels(xticklabels)
            ax[2 * p + 6 + d].set_yticklabels(yticklabels)
            cbar = heatmap.collections[0].colorbar
            cbar.ax.set_aspect(5)
            cbar.ax.yaxis.tick_left()
            cbar.ax.yaxis.set_label_position("left")
            if p == 0 and d == 0:
                ax[2 * p + 6 + d].set_ylabel("Fraction of trials")
            else:
                ax[2 * p + 6 + d].set_yticklabels([])
            if d == 0:
                ax[2 * p + 6 + d].set_title("$CSD$")
            else:
                ax[2 * p + 6 + d].set_title("$CSD^C$")
            if p == 1:
                x_min, x_max = ax[2 * p + 6 + d].get_xlim()
                y_min, y_max = ax[2 * p + 6 + d].get_ylim()
                rectangle = plt.Rectangle(
                    (1, 16),
                    x_max - 1.05,
                    y_max - 16,
                    facecolor="none",
                    edgecolor="r",
                    ls="--",
                    lw=2,
                )
                ax[2 * p + 6 + d].add_patch(rectangle)
                # ax[2 * p + 6 + d].plot(
                #     [1, len(bins) - 1], [16, 16], color="r", ls="--"
                # )
                # ax[2 * p + 6 + d].plot(
                #     [1, 1],
                #     [16, activation_data.shape[1]],
                #     color="r",
                #     ls="--",
                # )

    sns.despine()


def plot_legend(ax):
    """Plot legend

    Parameters
    ----------
    ax : axis object
    """
    ax.text(
        0.5,
        0.88,
        "Legend for B, C, and D, E, F, G",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    stimulus_driven_connected_trace_line = mlines.Line2D(
        [], [], ls="-", color="c", label="$CSD$"
    )
    stimulus_driven_not_connected_trace_line = mlines.Line2D(
        [], [], ls="-", color="m", label="$CSD^C$"
    )
    active_group_dot = mlines.Line2D(
        [], [], ls="", marker="o", color="k", alpha=0.5, label="$active$ $group$"
    )
    legend_handles = [
        stimulus_driven_connected_trace_line,
        stimulus_driven_not_connected_trace_line,
        active_group_dot,
    ]

    group_types = [
        "$CSD$",
        "$CSD^C$",
        "$CFM$",
        "$CFM^C$",
    ]
    colors = ["c", "m", "r", [0, 0, 0, 0.5]]
    for i, grp in enumerate(group_types):
        patch = mpatches.Patch(color=colors[i], label=grp)
        legend_handles.append(patch)
    ax.legend(handles=legend_handles, loc="center", frameon=False)
    ax.axis("off")


def main():
    figure = plt.figure(figsize=(13, 12), constrained_layout=True)
    spec = figure.add_gridspec(9, 15)
    axa = figure.add_subplot(spec[0:3, 3:9])
    # axb = figure.add_subplot(spec[0, 3:6])
    axc = figure.add_subplot(spec[0:3, 9:15])
    # axd = figure.add_subplot(spec[0, 9:12])
    axe = figure.add_subplot(spec[3:5, 3:6])
    axf = figure.add_subplot(spec[3:5, 6:9])
    axg = figure.add_subplot(spec[3:5, 9:12])
    axh = figure.add_subplot(spec[3:5, 12:15])
    axi = figure.add_subplot(spec[5:9, 3:6])
    axj = figure.add_subplot(spec[5:9, 6:9])
    axk = figure.add_subplot(spec[5:9, 9:12])
    axl = figure.add_subplot(spec[5:9, 12:15])
    axm = figure.add_subplot(spec[0:3, 0:3])
    axn = figure.add_subplot(spec[3:5, 0:3])
    axo = figure.add_subplot(spec[5:9, 0:3])

    plot_weibull(axm, M=4, a=10, powers=[2, 4, 8, 16])
    plot_activation_data(
        prm.hippo_elec,
        4,
        "weibull_exponent",
        [axa, axc, axe, axf, axg, axh, axi, axj, axk, axl, axo],
        powers=[2, 16],
    )
    plot_legend(axn)

    plt.figtext(0.02, 0.96, "A", fontsize=16, weight="bold")
    plt.figtext(0.24, 0.96, "B", fontsize=16, weight="bold")
    plt.figtext(0.63, 0.96, "C", fontsize=16, weight="bold")
    plt.figtext(0.23, 0.64, "D", fontsize=16, weight="bold")
    plt.figtext(0.42, 0.64, "E", fontsize=16, weight="bold")
    plt.figtext(0.61, 0.64, "F", fontsize=16, weight="bold")
    plt.figtext(0.80, 0.64, "G", fontsize=16, weight="bold")
    plt.figtext(0.23, 0.39, "H", fontsize=16, weight="bold")
    plt.figtext(0.42, 0.39, "I", fontsize=16, weight="bold")
    plt.figtext(0.61, 0.39, "J", fontsize=16, weight="bold")
    plt.figtext(0.80, 0.39, "K", fontsize=16, weight="bold")

    # Annotations
    plt.figtext(
        0.41,
        0.96,
        r"Weak nonlinearity, $\rho=2$",
        fontsize="large",
        horizontalalignment="center",
    )
    plt.figtext(
        0.81,
        0.96,
        r"Strong nonlinearity, $\rho=16$",
        fontsize="large",
        horizontalalignment="center",
    )
    figure.suptitle("Hippo-elec, 4 ensembles, group size=4")
    plt.savefig("Figure-3.png", bbox_inches="tight")
    plt.savefig("Figure-3.svg", bbox_inches="tight")

    plt.show()


if __name__ == "__main__":
    main()
