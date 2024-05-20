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
import seaborn as sns
import pandas as pd
from scipy.stats import mannwhitneyu


sys.path.append("../../codes/sequences/")
import sequence_params as prm

sys.path.append("../../../utils/")
from scalebars import add_scalebar

sns.set(style="ticks")
sns.set_context("paper")

DATA_PATH = "../../codes/sequences/activation_data/"

MIN_M = prm.MIN_M
MAX_M = 4

T_POP_POST = 4000
NUM_RUNS = 100


def read_patterns_data(
    net, dataset_name, trial_numbers, seq_len, num_patterns, datatype
):
    """Read sequence related data across trials

    Parameters
    ----------
    net : dict
        Dictionary of network parameters
    dataset_name : str
        Type of data file to be read
    trial_numbers : np.ndarray
        Array of trial numbers to read
    seq_len : int
        Sequence length
    num_patterns: int
        Number of patterns
    datatype:
        Datatype of the data to be read

    Returns
    -------
    data_matrix : np.ndarray
        3D array of data of shape num_trials x (t_pop_post x num_runs) x num_patterns
    """
    data_matrix = np.zeros(
        (len(trial_numbers), T_POP_POST * NUM_RUNS, num_patterns), dtype=datatype
    )
    for t, t_num in enumerate(trial_numbers):
        for r in range(NUM_RUNS):
            filename = DATA_PATH + "patterns-%s-%s-run-%03d-trial-%03d-M-%d.csv" % (
                dataset_name,
                net["name"],
                r + 1,
                t_num,
                seq_len,
            )
            data = np.asarray(pd.read_csv(filename, header=None))
            data_matrix[t, r * T_POP_POST : (r + 1) * T_POP_POST, :] = data

    return data_matrix


def convert_pvalue_to_stars(pvalue):
    """Convert pvalue to significance stars

    Parameters
    ----------
    pvalue : float
        p value obtained from the statistical tests

    Returns
    -------
        The right symbol for indicating the level of significance
    """
    if pvalue <= 0.0001:
        return "****"
    elif pvalue <= 0.001:
        return "***"
    elif pvalue <= 0.01:
        return "**"
    elif pvalue <= 0.05:
        return "*"
    return "ns"


def plot_activation_data(net, M, num_patterns, activation_type, ax):
    """Plot neuron activation data

    Parameters
    ----------
    net : dict
        Dictionary with network configuration parameters
    M : int
        Sequence length
    num_patterns : int
        Number of patterns
    activation : str
        String defining the type of activation
        "auc": Area under the curve
        "peak": Peak activation
    ax : list
        List of axis objects
    """
    stim_array = np.genfromtxt(
        DATA_PATH + f"patterns-stimulus-{net['name']}-run-001-M-{M}.csv", delimiter=","
    )

    connectivity_trial = np.where(stim_array < 0)[0] + 1
    cpsd_sequences = read_patterns_data(
        net=net,
        dataset_name="perfect-stimulus-driven-sequences",
        trial_numbers=connectivity_trial,
        seq_len=M,
        num_patterns=num_patterns,
        datatype=np.int32,
    )[0, :, 0]

    stim_trial_numbers = np.where(stim_array > 0)[0] + 1
    activation_data = read_patterns_data(
        net=net,
        dataset_name=f"activation-{activation_type}",
        trial_numbers=stim_trial_numbers,
        seq_len=M,
        num_patterns=num_patterns,
        datatype=np.float64,
    )
    pscsd_neurons_activation_data = activation_data[:, cpsd_sequences > 0, :]
    non_pscsd_neurons_activation_data = activation_data[:, cpsd_sequences == 0, :]

    bg_trial_numbers = np.where(stim_array == 0)[0] + 1
    bg_activation_data = read_patterns_data(
        net=net,
        dataset_name=f"activation-{activation_type}",
        trial_numbers=bg_trial_numbers,
        seq_len=M,
        num_patterns=1,
        datatype=np.float64,
    )
    pscsd_neurons_bg_activation_data = bg_activation_data[:, cpsd_sequences > 0, :]
    non_pscsd_neurons_bg_activation_data = bg_activation_data[:, cpsd_sequences == 0, :]

    all_sequences_stim_trials_data = read_patterns_data(
        net=net,
        dataset_name=f"all-sequences",
        trial_numbers=stim_trial_numbers,
        seq_len=M,
        num_patterns=num_patterns,
        datatype=np.int32,
    )

    patterns_df = pd.read_csv(
        f"{DATA_PATH}patterns-{net['name']}-run-001-M-{M}.csv", header=None
    )

    labels = ["$PSCSD$ $neurons$", "$PSCSD^C$ $neurons$"]
    min_actn = np.min(activation_data)
    max_actn = np.max(activation_data)
    gap = max_actn * 0.3

    for p, pattern_id in enumerate([0, int(num_patterns / 2), num_patterns - 1]):
        for n in np.arange(non_pscsd_neurons_activation_data.shape[1])[:3]:
            ax[p].plot(
                np.arange(non_pscsd_neurons_activation_data.shape[0]),
                n * gap + non_pscsd_neurons_activation_data[:, n, pattern_id],
                color="k",
                alpha=0.5,
                label="$PSCSD^C$ $neurons$",
            )
            hit_trials = np.where(
                all_sequences_stim_trials_data[:, cpsd_sequences == 0, pattern_id][:, n]
                > 0
            )[0]
            ax[p].scatter(
                hit_trials,
                n * gap + non_pscsd_neurons_activation_data[hit_trials, n, pattern_id],
                c="k",
                alpha=0.5,
                s=10,
                label="aPSD sequence",
            )

            if p == 0 and n == 0:
                ax[4].hist(
                    non_pscsd_neurons_activation_data[:, n, pattern_id],
                    bins=np.linspace(min_actn, max_actn, 11),
                    color="k",
                    alpha=0.5,
                    orientation="horizontal",
                )
                ax[4].set_ylabel("Activation")
                ax[4].set_xlabel("#Trials")
                ax[4].set_title("$PSCSD^C$ $Neuron$ $\#1$ $pattern$ $\#1$")

        for n in np.arange(pscsd_neurons_activation_data.shape[1])[:3]:
            ax[p].plot(
                np.arange(pscsd_neurons_activation_data.shape[0]),
                (4 + n) * gap + pscsd_neurons_activation_data[:, n, pattern_id],
                color="r",
                label="$PSCSD$ $neurons$",
            )
            hit_trials = np.where(
                all_sequences_stim_trials_data[:, cpsd_sequences > 0, pattern_id][:, n]
                > 0
            )[0]
            ax[p].scatter(
                hit_trials,
                (4 + n) * gap
                + pscsd_neurons_activation_data[hit_trials, n, pattern_id],
                c="k",
                alpha=0.5,
                s=10,
            )
            if p == 0 and n == 0:
                ax[3].hist(
                    pscsd_neurons_activation_data[:, n, pattern_id],
                    bins=np.linspace(min_actn, max_actn, 11),
                    color="r",
                    orientation="horizontal",
                )
                ax[3].set_ylabel("Activation")
                ax[3].set_xlabel("#Trials")
                ax[3].set_title("$PSCSD$ $Neuron$ $\#1$ $pattern$ $\#1$")
                neuron_types = [
                    "$PSCSD$ $neurons$",
                    "$PSCSD^C$ $neurons$",
                ]
                colors = ["r", [0, 0, 0, 0.5]]
                legend_handles = []
                for i, grp in enumerate(neuron_types):
                    patch = mpatches.Patch(color=colors[i], label=grp)
                    legend_handles.append(patch)
                ax[3].legend(
                    handles=legend_handles,
                    title="Legend for D,E,F,G,H,I,K",
                    loc="upper right",
                    frameon=False,
                )

        ax[p].set_axis_off()
        ax[p].set_title(
            f"Pattern #{pattern_id+1} - {list(patterns_df.iloc[pattern_id])}"
        )

        bins = np.linspace(
            np.min(np.mean(activation_data, axis=0)),
            np.max(np.mean(activation_data, axis=0)),
            11,
        )

        ax[5 + p].hist(
            [
                np.mean(pscsd_neurons_activation_data[:, :, pattern_id], axis=0),
                np.mean(non_pscsd_neurons_activation_data[:, :, pattern_id], axis=0),
            ],
            bins=bins,
            log=True,
            color=["r", "#00000080"],
            label=labels,
        )
        ax[5 + p].set_xlabel("Trial averaged activation")
        if p == 0:
            ax[5 + p].set_ylabel("#Neurons")
        else:
            ax[5 + p].set_yticklabels([])

        if p == 0:
            perfect_stimulus_driven_connected_trace_line = mlines.Line2D(
                [], [], ls="-", color="r", label="$PSCSD$"
            )
            perfect_stimulus_driven_not_connected_trace_line = mlines.Line2D(
                [], [], ls="-", color="k", alpha=0.5, label="$PSCSD^C$"
            )
            # perfect_active_sequence_dot = mlines.Line2D(
            #     [], [], ls="", marker="o", color="r", label="$active$ $PSD$ $sequence$"
            # )
            active_sequence_dot = mlines.Line2D(
                [],
                [],
                ls="",
                marker="o",
                color="k",
                alpha=0.5,
                label="$active$ $sequence$",
            )
            legend_handles = [
                perfect_stimulus_driven_connected_trace_line,
                perfect_stimulus_driven_not_connected_trace_line,
                active_sequence_dot,
            ]
            ax[p].legend(handles=legend_handles, loc="center", frameon=False)
            ax[p].axis("off")

        if p == 1:
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

    pscsd_selectivities = np.zeros(
        pscsd_neurons_activation_data.shape[1], dtype=np.float32
    )
    for neuron in range(20):
        ax[9].plot(
            np.arange(1, num_patterns + 1),
            np.mean(pscsd_neurons_activation_data[:, neuron, :], axis=0),
        )
        pscsd_selectivities[neuron] = (
            np.mean(pscsd_neurons_activation_data[:, neuron, 0], axis=0)
            - np.mean(np.mean(pscsd_neurons_activation_data[:, neuron, :], axis=0))
        ) / np.max(np.mean(pscsd_neurons_activation_data[:, neuron, :], axis=0))

    ax[9].set_xlabel("Pattern #")
    ax[9].set_ylabel("Trial averaged activation")
    ax[9].set_title("PSCSD neurons")

    non_pscsd_selectivities = np.zeros(
        non_pscsd_neurons_activation_data.shape[1], dtype=np.float32
    )
    for neuron in range(non_pscsd_neurons_activation_data.shape[1]):
        non_pscsd_selectivities[neuron] = (
            np.mean(non_pscsd_neurons_activation_data[:, neuron, 0], axis=0)
            - np.mean(np.mean(non_pscsd_neurons_activation_data[:, neuron, :], axis=0))
        ) / np.max(np.mean(non_pscsd_neurons_activation_data[:, neuron, :], axis=0))

    print(
        f"PSCSD, min={min(pscsd_selectivities)}, max={max(pscsd_selectivities)}, mean={np.mean(pscsd_selectivities)}, std={np.std(pscsd_selectivities)}"
    )

    print(
        f"non_pscsd, min={min(non_pscsd_selectivities)}, max={max(non_pscsd_selectivities)}, mean={np.mean(non_pscsd_selectivities)}, std={np.std(non_pscsd_selectivities)}"
    )
    _, p_val = mannwhitneyu(
        pscsd_selectivities, non_pscsd_selectivities, alternative="greater"
    )
    print(f"{p_val=}")

    pscsd_bg_subtracted_activation_trial_averaged = np.mean(
        pscsd_neurons_activation_data, axis=0
    ) - np.mean(pscsd_neurons_bg_activation_data, axis=0)
    non_pscsd_bg_subtracted_activation_trial_averaged = np.mean(
        non_pscsd_neurons_activation_data, axis=0
    ) - np.mean(non_pscsd_neurons_bg_activation_data, axis=0)

    pscsd_selectivities_bg_subtracted = np.zeros(
        pscsd_bg_subtracted_activation_trial_averaged.shape[0], dtype=np.float32
    )
    non_pscsd_selectivities_bg_subtracted = np.zeros(
        non_pscsd_bg_subtracted_activation_trial_averaged.shape[0], dtype=np.float32
    )
    bins = np.linspace(
        np.min(np.mean(bg_activation_data, axis=0)),
        np.max(np.mean(bg_activation_data, axis=0)),
        11,
    )

    ax[8].hist(
        [
            np.mean(pscsd_neurons_bg_activation_data[:, :, 0], axis=0),
            np.mean(non_pscsd_neurons_bg_activation_data[:, :, 0], axis=0),
        ],
        bins=bins,
        log=True,
        color=["r", "#00000080"],
    )
    ax[8].set_xlabel("Trial averaged activation")
    ax[8].set_ylabel("#Neurons")
    ax[8].set_title("Background activity trials")
    for neuron in range(pscsd_bg_subtracted_activation_trial_averaged.shape[0]):
        pscsd_selectivities_bg_subtracted[neuron] = (
            pscsd_bg_subtracted_activation_trial_averaged[neuron, 0]
            - np.mean(pscsd_bg_subtracted_activation_trial_averaged[neuron, :])
        ) / np.max(pscsd_bg_subtracted_activation_trial_averaged[neuron, :])

    for neuron in range(non_pscsd_bg_subtracted_activation_trial_averaged.shape[0]):
        non_pscsd_selectivities_bg_subtracted[neuron] = (
            non_pscsd_bg_subtracted_activation_trial_averaged[neuron, 0]
            - np.mean(non_pscsd_bg_subtracted_activation_trial_averaged[neuron, :])
        ) / np.max(non_pscsd_bg_subtracted_activation_trial_averaged[neuron, :])

    print(
        f"PSCSD, min={min(pscsd_selectivities_bg_subtracted)}, max={max(pscsd_selectivities_bg_subtracted)}, mean={np.mean(pscsd_selectivities_bg_subtracted)}, std={np.std(pscsd_selectivities_bg_subtracted)}"
    )
    print(
        f"non_pscsd, min={min(non_pscsd_selectivities_bg_subtracted)}, max={max(non_pscsd_selectivities_bg_subtracted)}, mean={np.mean(non_pscsd_selectivities_bg_subtracted)}, std={np.std(non_pscsd_selectivities_bg_subtracted)}"
    )
    _, p_val_bg_subtracted = mannwhitneyu(
        pscsd_selectivities_bg_subtracted,
        non_pscsd_selectivities_bg_subtracted,
        alternative="greater",
    )
    print(f"{p_val_bg_subtracted=}")

    violins = ax[10].violinplot(
        [
            pscsd_selectivities,
            non_pscsd_selectivities,
            pscsd_selectivities_bg_subtracted,
            non_pscsd_selectivities_bg_subtracted,
        ],
        showmedians=True,
    )
    violin_colors = ["r", "k", "r", "k"]
    alphas = [0.5, 0.25, 0.5, 0.25]
    for i, pc in enumerate(violins["bodies"]):
        pc.set_facecolor(violin_colors[i])
        pc.set_alpha(alphas[i])
    violins["cbars"].set_colors(violin_colors)
    violins["cmins"].set_colors(violin_colors)
    violins["cmaxes"].set_colors(violin_colors)
    violins["cmedians"].set_colors(violin_colors)
    ax[10].set_xticks(
        np.arange(1, len(labels) * 2 + 1),
        [
            "$PSCSD$",
            "$PSCSD^C$",
            "$PSCSD$-$bg$",
            "$PSCSD^C$-$bg$",
        ],
        rotation=45,
    )
    y_max = max(np.max(pscsd_selectivities), np.max(non_pscsd_selectivities))
    ax[10].plot(
        [1, 2],
        [y_max + 0.015, y_max + 0.015],
        color="k",
    )
    ax[10].text(
        x=1.5,
        y=y_max + 0.03,
        s=convert_pvalue_to_stars(p_val),
        va="center",
        ha="center",
        fontsize=8,
    )
    y_max = max(
        np.max(pscsd_selectivities_bg_subtracted),
        np.max(non_pscsd_selectivities_bg_subtracted),
    )
    ax[10].plot(
        [3, 4],
        [y_max + 0.015, y_max + 0.015],
        color="k",
    )
    ax[10].text(
        x=3.5,
        y=y_max + 0.03,
        s=convert_pvalue_to_stars(p_val_bg_subtracted),
        va="center",
        ha="center",
        fontsize=8,
    )
    ax[10].set_ylabel("Selectivity")

    sns.despine()


def main():
    figure = plt.figure(figsize=(12, 13), constrained_layout=True)
    spec = figure.add_gridspec(12, 16)
    axa = figure.add_subplot(spec[0:6, 0:4])
    axb = figure.add_subplot(spec[0:6, 4:8])
    axc = figure.add_subplot(spec[0:6, 8:12])
    axd = figure.add_subplot(spec[0:3, 12:16])
    axe = figure.add_subplot(spec[3:6, 12:16])
    axf = figure.add_subplot(spec[6:9, 0:4])
    axg = figure.add_subplot(spec[6:9, 4:8])
    axh = figure.add_subplot(spec[6:9, 8:12])
    axi = figure.add_subplot(spec[6:9, 12:16])
    axj = figure.add_subplot(spec[9:12, 0:12])
    axk = figure.add_subplot(spec[9:12, 12:16])

    plot_activation_data(
        net=prm.cortex_cicr,
        M=4,
        num_patterns=24,
        activation_type="auc",
        ax=[
            axa,
            axb,
            axc,
            axd,
            axe,
            axf,
            axg,
            axh,
            axi,
            axj,
            axk,
        ],
    )

    plt.figtext(0.04, 0.97, "A", fontsize=16, weight="bold")
    plt.figtext(0.27, 0.97, "B", fontsize=16, weight="bold")
    plt.figtext(0.51, 0.97, "C", fontsize=16, weight="bold")
    plt.figtext(0.75, 0.97, "D", fontsize=16, weight="bold")
    plt.figtext(0.75, 0.71, "E", fontsize=16, weight="bold")
    plt.figtext(0.02, 0.46, "F", fontsize=16, weight="bold")
    plt.figtext(0.25, 0.46, "G", fontsize=16, weight="bold")
    plt.figtext(0.48, 0.46, "H", fontsize=16, weight="bold")
    plt.figtext(0.73, 0.46, "I", fontsize=16, weight="bold")
    plt.figtext(0.01, 0.22, "J", fontsize=16, weight="bold")
    plt.figtext(0.73, 0.22, "K", fontsize=16, weight="bold")

    figure.suptitle("Cortex-CICR, 4 ensembles, sequence length=4")
    plt.savefig("Figure-7-cortex-cicr.png", bbox_inches="tight")
    plt.savefig("Figure-7-cortex-cicr.svg", bbox_inches="tight")

    plt.show()


if __name__ == "__main__":
    main()
