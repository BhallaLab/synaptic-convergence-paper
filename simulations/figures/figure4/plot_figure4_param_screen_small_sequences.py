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


def plot_prob_for_a_fraction(ax, net, M):
    """Plot parameter scan for probability of c_poss

    Parameters
    ----------
    ax : axis object
    net : dict
        Dictionary with network configuration parameters
    M : int
        Sequence length, also equal to the number of ensembles
    """
    p_values = np.logspace(-3, 0, num=7)
    N_values = np.logspace(0, 4, num=9)
    L = net["dendrite_length"]
    delta = net["delta"]

    # connectivity-based perfectly ordered stimulus sequences
    p_cposs_seq = np.zeros((len(N_values), len(p_values)))
    for n, N in enumerate(N_values):
        for p_i, p in enumerate(p_values):
            E_cposs = p * N * np.power(p * N * delta / L, M - 1)
            p_cposs_seq[n, p_i] = 1 - np.exp(-E_cposs)

    cf = ax.contour(
        p_values,
        N_values,
        p_cposs_seq,
        levels=np.logspace(-8, 0, 9),
        cmap=cm.viridis,
        norm=LogNorm(1e-8, 1e0),
    )
    cbar = plt.colorbar(
        cf,
        ax=ax,
        format="%.0e",
    )
    cbar.ax.set_title("p_cposs")
    print(np.min(p_cposs_seq))
    print(np.max(p_cposs_seq))
    cbar.set_ticks([1e-8, 1e-6, 1e-4, 1e-2, 1e0])
    sns.despine()
    ax.set_xlabel("Connection probability p")
    ax.set_ylabel("Ensemble size N")
    ax.set_xscale("log", nonpositive="mask")
    ax.set_yscale("log", nonpositive="mask")


def main():
    figure, ax = plt.subplots(1, 2, figsize=(9, 4), constrained_layout=True)

    plot_prob_for_a_fraction(ax[0], networks[0], 2)
    ax[0].set_title("Chem & CICR")
    plot_prob_for_a_fraction(ax[1], networks[2], 2)
    ax[1].set_title("Elec")

    plt.savefig("Figure-4-param-screen-small-sequences.png", bbox_inches="tight")
    plt.savefig("Figure-4-param-screen-small-sequences.svg", bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
