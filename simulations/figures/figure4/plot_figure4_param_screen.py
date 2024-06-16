#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : plot_figure.py
# Author            : Bhanu Priya S <bhanups@ncbs.res.in>
# Date              : 01.02.2023
# Last Modified Date: 11.05.2023
# Last Modified By  : Bhanu Priya S <bhanups@ncbs.res.in>

import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.ticker as mticker
import seaborn as sns
from matplotlib import cm
from matplotlib.colors import LogNorm


sys.path.append("../../codes/sequences/")
import sequence_params as prm

sns.set(style="ticks")
sns.set_context("paper")


networks = [
    prm.hippo_chem,
    prm.hippo_cicr,
    prm.hippo_elec,
    prm.cortex_chem,
    prm.cortex_cicr,
    prm.cortex_elec,
]
colors = sns.color_palette("muted", n_colors=len(networks))


def plot_conn_prob_vs_p_seq(ax, net, M):
    """Plot probability of different sequences as a function of connection probability

    Parameters
    ----------
    ax : axis object
    net : dict
        Dictionary with network configuration parameters
    M : int
        Sequence length, also equal to the number of ensembles
    """
    # Parameters of cortex_cicr
    p = np.logspace(-2, 0, num=7)
    p_e = net["ensemble_participation_prob"]
    N = net["ensemble_size"]
    L = net["dendrite_length"]
    R = net["background_rate"]
    D = net["sequence_time_step"]
    t_pop_pre = net["population_size"]
    delta = net["delta"]
    bg_prob = 1 - np.exp(-R * D)
    num_synapses = np.array([int(j * t_pop_pre) for j in p])

    E_cposs = p * N * np.power(p * N * delta / L, M - 1)
    E_aposs = p * p_e * N * np.power(p * p_e * N * delta / L, M - 1)
    E_all_seq = (p * p_e * N + bg_prob * (num_synapses - p * N)) * np.power(
        (p * p_e * N + bg_prob * (num_synapses - p * N)) * delta / L, M - 1
    )
    E_noise_seq = (
        bg_prob
        * (num_synapses - p * N)
        * np.power(bg_prob * (num_synapses - p * N) * delta / L, M - 1)
    )
    E_gap_fill_seq = E_all_seq - E_noise_seq - E_aposs

    p_cposs_seq = 1 - np.exp(-E_cposs)
    p_aposs_seq = 1 - np.exp(-E_aposs)
    p_noise_seq = 1 - np.exp(-E_noise_seq)
    p_gap_fill_seq = 1 - np.exp(-E_gap_fill_seq)
    p_any_seq = 1 - np.exp(-E_all_seq)

    ax.plot(p, p_cposs_seq, "-o", c="b", label="cPOSS")
    ax.plot(p, p_aposs_seq, ":o", c="b", label="aPOSS")
    ax.plot(p, p_gap_fill_seq, ":s", c="c", label="GS")
    ax.plot(p, p_noise_seq, ":+", c="r", label="NS")
    ax.plot(p, p_any_seq, ":x", c="k", label="AS")

    sns.despine()
    ax.set_xscale("log", nonpositive="mask")
    ax.set_yscale("log", nonpositive="mask")
    ax.set_ylim([1e-5, 1e0])
    ax.set_xlabel("Connection probability, p")
    ax.set_ylabel("p_seq")
    ax.set_xticks([1e-2, 1e-1, 1e0])
    ax.legend(title="Legend for A,B,C,D", frameon=False)


def plot_ensemble_participation_prob_vs_p_group(ax, net, M):
    """Plot probability of different sequences as a function of ensemble participation probability

    Parameters
    ----------
    ax : axis object
    net : dict
        Dictionary with network configuration parameters
    M : int
        Sequence length, also equal to the number of ensembles
    """
    # Parameters of cortex_cicr
    p = net["connection_prob"]
    p_e = np.linspace(0, 1, 11)
    N = net["ensemble_size"]
    L = net["dendrite_length"]
    R = net["background_rate"]
    D = net["sequence_time_step"]
    t_pop_pre = net["population_size"]
    delta = net["delta"]
    bg_prob = 1 - np.exp(-R * D)
    num_synapses = int(p * t_pop_pre)

    E_cposs = p * N * np.power(p * N * delta / L, M - 1)
    E_aposs = p * p_e * N * np.power(p * p_e * N * delta / L, M - 1)
    E_all_seq = (p * p_e * N + bg_prob * (num_synapses - p * N)) * np.power(
        (p * p_e * N + bg_prob * (num_synapses - p * N)) * delta / L, M - 1
    )
    E_noise_seq = (
        bg_prob
        * (num_synapses - p * N)
        * np.power(bg_prob * (num_synapses - p * N) * delta / L, M - 1)
    )
    E_gap_fill_seq = E_all_seq - E_noise_seq - E_aposs

    p_cposs_seq = [1 - np.exp(-E_cposs)] * len(p_e)
    p_aposs_seq = 1 - np.exp(-E_aposs)
    p_noise_seq = [1 - np.exp(-E_noise_seq)] * len(p_e)
    p_gap_fill_seq = 1 - np.exp(-E_gap_fill_seq)
    p_any_seq = 1 - np.exp(-E_all_seq)

    ax.plot(p_e, p_cposs_seq, "-o", c="b", label="cPOSS")
    ax.plot(p_e, p_aposs_seq, ":o", c="b", label="aPOSS")
    ax.plot(p_e, p_gap_fill_seq, ":s", c="c", label="GS")
    ax.plot(p_e, p_noise_seq, ":+", c="r", label="NS")
    ax.plot(p_e, p_any_seq, ":x", c="k", label="AS")

    sns.despine()
    ax.set_yscale("log", nonpositive="mask")
    ax.set_ylim([1e-5, 1e0])
    ax.set_xlabel("Ensemble participation\nprobability, p_e")
    ax.set_ylabel("p_seq")


def plot_input_zone_width_vs_p_seq(ax, net, M):
    """Plot probability of different sequences as a function of input zone width

    Parameters
    ----------
    ax : axis object
    net : dict
        Dictionary with network configuration parameters
    M : int
        Sequence length, also equal to the number of ensembles
    """
    # Parameters of cortex_cicr
    p = net["connection_prob"]
    p_e = net["ensemble_participation_prob"]
    N = net["ensemble_size"]
    L = net["dendrite_length"]
    R = net["background_rate"]
    D = net["sequence_time_step"]
    t_pop_pre = net["population_size"]
    delta = np.linspace(1, 20, 20) * 1e-6
    bg_prob = 1 - np.exp(-R * D)
    num_synapses = int(p * t_pop_pre)

    E_cposs = p * N * np.power(p * N * delta / L, M - 1)
    E_aposs = p * p_e * N * np.power(p * p_e * N * delta / L, M - 1)
    E_all_seq = (p * p_e * N + bg_prob * (num_synapses - p * N)) * np.power(
        (p * p_e * N + bg_prob * (num_synapses - p * N)) * delta / L, M - 1
    )
    E_noise_seq = (
        bg_prob
        * (num_synapses - p * N)
        * np.power(bg_prob * (num_synapses - p * N) * delta / L, M - 1)
    )
    E_gap_fill_seq = E_all_seq - E_noise_seq - E_aposs

    p_cposs_seq = 1 - np.exp(-E_cposs)
    p_aposs_seq = 1 - np.exp(-E_aposs)
    p_noise_seq = 1 - np.exp(-E_noise_seq)
    p_gap_fill_seq = 1 - np.exp(-E_gap_fill_seq)
    p_any_seq = 1 - np.exp(-E_all_seq)

    ax.plot(delta, p_cposs_seq, "-o", c="b", label="cPOSS")
    ax.plot(delta, p_aposs_seq, ":o", c="b", label="aPOSS")
    ax.plot(delta, p_gap_fill_seq, ":s", c="c", label="GS")
    ax.plot(delta, p_noise_seq, ":+", c="r", label="NS")
    ax.plot(delta, p_any_seq, ":x", c="k", label="AS")

    sns.despine()
    ax.set_yscale("log", nonpositive="mask")
    ax.set_ylim([1e-5, 1e0])
    ax.set_xlabel(r"Input zone width, $\Delta$")
    ax.set_ylabel("p_seq")


def plot_ensemble_size_vs_p_seq(ax, net, M):
    """Plot probability of different sequences as a function of ensemble size

    Parameters
    ----------
    ax : axis object
    net : dict
        Dictionary with network configuration parameters
    M : int
        Sequence length, also equal to the number of ensembles
    """
    # Parameters of cortex_cicr
    p = net["connection_prob"]
    p_e = net["ensemble_participation_prob"]
    N = np.logspace(1, 4, num=7)
    L = net["dendrite_length"]
    R = net["background_rate"]
    D = net["sequence_time_step"]
    t_pop_pre = net["population_size"]
    delta = net["delta"]
    bg_prob = 1 - np.exp(-R * D)
    num_synapses = int(p * t_pop_pre)

    E_cposs = p * N * np.power(p * N * delta / L, M - 1)
    E_aposs = p * p_e * N * np.power(p * p_e * N * delta / L, M - 1)
    E_all_seq = (p * p_e * N + bg_prob * (num_synapses - p * N)) * np.power(
        (p * p_e * N + bg_prob * (num_synapses - p * N)) * delta / L, M - 1
    )
    E_noise_seq = (
        bg_prob
        * (num_synapses - p * N)
        * np.power(bg_prob * (num_synapses - p * N) * delta / L, M - 1)
    )
    E_gap_fill_seq = E_all_seq - E_noise_seq - E_aposs

    p_cposs_seq = 1 - np.exp(-E_cposs)
    p_aposs_seq = 1 - np.exp(-E_aposs)
    p_noise_seq = 1 - np.exp(-E_noise_seq)
    p_gap_fill_seq = 1 - np.exp(-E_gap_fill_seq)
    p_any_seq = 1 - np.exp(-E_all_seq)

    ax.plot(N, p_cposs_seq, "-o", c="b", label="cPOSS")
    ax.plot(N, p_aposs_seq, ":o", c="b", label="aPOSS")
    ax.plot(N, p_gap_fill_seq, ":s", c="c", label="GS")
    ax.plot(N, p_noise_seq, ":+", c="r", label="NS")
    ax.plot(N, p_any_seq, ":x", c="k", label="AS")

    sns.despine()
    ax.set_xscale("log", nonpositive="mask")
    ax.set_yscale("log", nonpositive="mask")
    ax.set_ylim([1e-5, 1e0])
    ax.set_xlabel("Ensemble size, N")
    ax.set_xticks([1e1, 1e2, 1e3, 1e4])
    ax.set_ylabel("p_seq")


def plot_p_noise_seq_vs_noise(ax, net, M):
    """Plot probability of any sequence as a function of the rate of background activity and input duration

    Parameters
    ----------
    ax : axis object
    net : dict
        Dictionary with network configuration parameters
    M : int
        Sequence length, also equal to the number of ensembles
    """
    # Parameters of cortex_cicr
    p = net["connection_prob"]
    N = net["ensemble_size"]
    L = net["dendrite_length"]
    R = np.logspace(-3, 1, num=9)
    D = np.logspace(-3, 1, num=9)
    t_pop_pre = net["population_size"]
    delta = net["delta"]
    X, Y = np.meshgrid(R, D)
    bg_prob = 1 - np.exp(-X * Y)
    num_synapses = int(p * t_pop_pre)

    E_noise_seq = (
        bg_prob
        * (num_synapses - p * N)
        * np.power(bg_prob * (num_synapses - p * N) * delta / L, M - 1)
    )

    p_noise_seq = 1 - np.exp(-E_noise_seq)
    p_noise_seq = np.clip(p_noise_seq, a_min=1e-6, a_max=None)

    ax = sns.heatmap(
        p_noise_seq,
        ax=ax,
        norm=LogNorm(p_noise_seq.min(), p_noise_seq.max()),
        xticklabels=R,
        yticklabels=D,
        cmap=cm.viridis,
        cbar_kws={
            "extend": "min",
            "ticks": [1e-6, 1e-3, 1e0],
            "format": mticker.FixedFormatter(["<=$10^{-6}$", "$10^{-3}$", "$10^{0}$"]),
        },
    )
    cbar = ax.collections[0].colorbar
    cbar.ax.set_title("p_noise_seq")

    fmt = "{:0.3f}"
    xticklabels = []
    for item in ax.get_xticklabels():
        item.set_text(fmt.format(float(item.get_text())))
        xticklabels += [item]

    yticklabels = []
    for item in ax.get_yticklabels():
        item.set_text(fmt.format(float(item.get_text())))
        yticklabels += [item]

    sns.despine()
    ax.set_xticklabels(xticklabels)
    ax.set_yticklabels(yticklabels)
    ax.invert_yaxis()
    ax.set_xlabel("Background activity(Hz)")
    ax.set_ylabel("Input duration(s)")


def plot_p_gapfill_vs_p_fp_seq(ax, net, M):
    """Plot fraction of gap-fill sequences amongst all false positive sequences for different ensemble sizes

    Parameters
    ----------
    ax : axis object
    net : dict
        Dictionary with network configuration parameters
    M : int
        Sequence length, also equal to the number of ensembles
    """
    scaling_factor = 8

    p = net["connection_prob"]
    p_e = net["ensemble_participation_prob"]
    N = np.logspace(1, 4, num=7)
    L = net["dendrite_length"]
    M = 5
    R = net["background_rate"]
    D = net["sequence_time_step"]
    t_pop_pre = net["population_size"]
    delta = net["delta"]
    bg_prob = 1 - np.exp(-R * D)
    num_synapses = int(p * t_pop_pre)

    E_aposs = p * p_e * N * np.power(p * p_e * N * delta / L, M - 1)
    E_all_seq = (p * p_e * N + bg_prob * (num_synapses - p * N)) * np.power(
        (p * p_e * N + bg_prob * (num_synapses - p * N)) * delta / L, M - 1
    )
    E_noise_seq = (
        bg_prob
        * (num_synapses - p * N)
        * np.power(bg_prob * (num_synapses - p * N) * delta / L, M - 1)
    )
    E_gap_fill_seq = E_all_seq - E_noise_seq - E_aposs

    E_fp_seq = E_all_seq - E_aposs

    ax.scatter(E_fp_seq, E_gap_fill_seq, c="r", s=np.log(N) * scaling_factor)
    ax.plot(E_fp_seq, E_gap_fill_seq, c="r", ls="--")

    sns.despine(ax=ax)
    x = np.logspace(-8, 4, num=7)
    ax.plot(x, x, color="k", ls="-", alpha=0.3)
    ax.set_xscale("log", nonpositive="mask")
    ax.set_yscale("log", nonpositive="mask")
    ax.set_xlim([1e-8, 1e4])
    ax.set_ylim([1e-8, 1e4])
    ax.set_xlabel("E_false_positive_seq")
    ax.set_ylabel("E_gapfill_seq")
    legend_handles = []
    for e_size in N:
        point = mlines.Line2D(
            [],
            [],
            marker="o",
            markersize=np.sqrt(np.log(e_size) * scaling_factor),
            ls="",
            color="r",
            label=f"{round(e_size)}",
        )
        legend_handles.append(point)
    ax.legend(handles=legend_handles, frameon=False)


def main():
    figure = plt.figure(figsize=(7, 8), constrained_layout=True)
    spec = figure.add_gridspec(3, 2)
    axa = figure.add_subplot(spec[0, 0])
    axb = figure.add_subplot(spec[0, 1])
    axc = figure.add_subplot(spec[1, 0])
    axd = figure.add_subplot(spec[1, 1])
    axe = figure.add_subplot(spec[2, 0])
    axf = figure.add_subplot(spec[2, 1])

    plot_conn_prob_vs_p_seq(axa, networks[4], 5)
    plot_ensemble_participation_prob_vs_p_group(axb, networks[4], 5)
    plot_input_zone_width_vs_p_seq(axc, networks[4], 5)
    plot_ensemble_size_vs_p_seq(axd, networks[4], 5)
    plot_p_noise_seq_vs_noise(axe, networks[4], 5)
    plot_p_gapfill_vs_p_fp_seq(axf, networks[4], 5)

    plt.figtext(0.01, 0.98, "A", fontsize=12, weight="bold")
    plt.figtext(0.55, 0.98, "B", fontsize=12, weight="bold")
    plt.figtext(0.01, 0.66, "C", fontsize=12, weight="bold")
    plt.figtext(0.55, 0.66, "D", fontsize=12, weight="bold")
    plt.figtext(0.01, 0.32, "E", fontsize=12, weight="bold")
    plt.figtext(0.55, 0.32, "F", fontsize=12, weight="bold")

    plt.savefig("Figure-4-param-screen.png")
    plt.savefig("Figure-4-param-screen.svg")

    plt.show()


if __name__ == "__main__":
    main()
