#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : plot_figure2_param_screen.py
# Author            : Bhanu Priya S <bhanups@ncbs.res.in>
# Date              : 09.11.2023
# Last Modified Date: 09.11.2023
# Last Modified By  : Bhanu Priya S <bhanups@ncbs.res.in>

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm
import seaborn as sns


sys.path.append("../../codes/groups/")
import group_params as prm

sns.set(style="ticks")
sns.set_context("paper")

DATA_PATH = "../../codes/groups/data/"

networks = [
    prm.hippo_chem,
    prm.hippo_cicr,
    prm.hippo_elec,
    prm.cortex_chem,
    prm.cortex_cicr,
    prm.cortex_elec,
]


def plot_prob_for_a_fraction(ax, net, m):
    """Plot parameter scan for probability of cFMG

    Parameters
    ----------
    ax : axis object
    net : dict
        Dictionary with network configuration parameters
    m : int
        Group size
    """
    p_values = np.logspace(-3, 0, num=7)
    N_values = np.logspace(0, 4, num=9)
    L = net["dendrite_length"]
    Z = net["zone_length"]

    # connectivity-based fully-mixed group
    c = 2.56141683
    p_total_exp_fit = np.zeros((len(N_values), len(p_values)))
    for n, N in enumerate(N_values):
        for p_i, p in enumerate(p_values):
            egrp_exp_Ei = p * N * Z / L
            p_grp_exp = np.power(1 - np.exp(-egrp_exp_Ei), m)
            p_total_exp_fit[n, p_i] = 1 - np.power((1 - p_grp_exp), c * int(L / Z))

    cf = ax.contour(
        p_values,
        N_values,
        p_total_exp_fit,
        levels=np.logspace(-8, 0, 9),
        cmap=cm.viridis,
        norm=LogNorm(1e-8, 1e0),
    )
    cbar = plt.colorbar(
        cf,
        ax=ax,
        format="%.0e",
    )
    cbar.ax.set_title("p_cFMG")
    print(np.min(p_total_exp_fit))
    print(np.max(p_total_exp_fit))
    cbar.set_ticks([1e-8, 1e-6, 1e-4, 1e-2, 1e0])
    sns.despine()
    ax.set_xlabel("Connection probability p")
    ax.set_ylabel("Ensemble size N")
    ax.set_xscale("log", nonpositive="mask")
    ax.set_yscale("log", nonpositive="mask")


def main():
    figure, ax = plt.subplots(1, 2, figsize=(9, 4), constrained_layout=True)
    plt.figtext(0.01, 0.96, "A", fontsize=12, weight="bold")
    plt.figtext(0.50, 0.96, "B", fontsize=12, weight="bold")

    plot_prob_for_a_fraction(ax[0], networks[0], 2)
    ax[0].set_title("Chem & CICR")
    plot_prob_for_a_fraction(ax[1], networks[2], 2)
    ax[1].set_title("Elec")

    plt.savefig("Figure-2-param-screen-small-sequences.png", bbox_inches="tight")
    plt.savefig("Figure-2-param-screen-small-sequences.svg", bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
