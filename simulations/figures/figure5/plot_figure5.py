#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : plot_figure6.py
# Author            : Bhanu Priya S <bhanups@ncbs.res.in>
# Date              : 01.02.2023
# Last Modified Date: 09.11.2023
# Last Modified By  : Bhanu Priya S <bhanups@ncbs.res.in>

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import seaborn as sns
import pandas as pd
from scipy.stats import wilcoxon

sns.set(style="ticks")
sns.set_context("paper")

DATA_PATH = "../../codes/sequences/chemical/data/"


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


def plot_panel_a(ax):
    """Display schematic of bistable-switch model

    Parameters
    ----------
    ax : axis object
    """
    ax.imshow(mpimg.imread("bistable_switch.png"))
    ax.axis("off")


def plot_panel_b(ax):
    # TODO put actual model
    """Plot response of bistable-switch model to sequential vs scrambled inputs

    Parameters
    ----------
    ax : axis object
    """
    ordered_response = pd.read_csv(
        DATA_PATH + "ordered_output.csv",
        header=None,
        sep=",",
    )
    scrambled_response = pd.read_csv(
        DATA_PATH + "scrambled_output.csv",
        header=None,
        sep=",",
    )
    start_ind = 25
    end_ind = 75
    ax.plot(
        np.arange(start_ind, end_ind),
        ordered_response.iloc[0][start_ind:end_ind],
        label="sequence",
    )
    ax.plot(
        np.arange(start_ind, end_ind),
        scrambled_response.iloc[0][start_ind:end_ind],
        label="scrambled",
    )
    ax.legend(frameon=False, loc="upper left")
    ax.set_xlabel("Position (\u03bcm)")
    ax.set_ylabel("A (a.u.)")
    sns.despine()


def plot_panel_c(ax):
    """Plot responses to different stimulus patterns as a function of the directionality of stimulus

    Parameters
    ----------
    ax : axis object
    """
    data = np.load(
        DATA_PATH + "responses_to_stim_patterns.npy",
        allow_pickle=True,
    ).item()
    ax.scatter(
        data["Q_score"][1:-1],
        data["A_total"][1:-1],
        label="other patterns",
        color="C1",
    )
    ax.scatter(
        data["Q_score"][0],
        data["A_total"][0],
        marker="x",
        color="C0",
        label="ordered sequence",
    )
    ax.scatter(
        data["Q_score"][-1],
        data["A_total"][-1],
        marker="^",
        color="C0",
        label="reverse sequence",
    )
    p = np.poly1d(np.polyfit(data["Q_score"], data["A_total"], 6))
    x_points = np.linspace(-1, 1, 30)
    ax.plot(x_points, p(x_points), color="k", ls="--", label="polynomial fit")
    ax.set_xlabel("Directionality (Q score)")
    ax.set_ylabel("A total (a.u.)")
    ax.set_xlim([-1.05, 1.05])
    ax.set_ylim(ymin=0)
    ax.legend(frameon=False)
    sns.despine()


def plot_panel_d(ax):
    ax.imshow(mpimg.imread("ectopic_inputs.png"))
    ax.axis("off")


def plot_panel_e(ax):
    """Plot selectivity in the absence and presence of ectopic inputs

    Parameters
    ----------
    ax : axis object
    """
    data = np.load(
        DATA_PATH + "ectopic_positions_selectivity.npy",
        allow_pickle=True,
    ).item()
    ax.axhline(data["reference"], color="k", label="reference")
    ax.plot(data["positions"], data["t_start"], "-", label="t_start")
    ax.plot(data["positions"], data["t_mid"], "-", label="t_mid")
    ax.plot(data["positions"], data["t_end"], "-", label="t_end")
    ax.axvspan(data["zone"][0], data["zone"][-1], color="c", alpha=0.5, lw=0)
    ax.legend(frameon=False)
    ax.set_xlabel("Position of extra input (\u03bcm)")
    ax.set_ylabel("Selectivity")
    sns.despine()


def plot_panel_f(ax):
    """Plot selectivity in the absence of ectopic inputs, and in the presence of 1 or 2 ectopic inputs

    Parameters
    ----------
    ax : axis object
    """
    data = np.load(
        DATA_PATH + "seq_lengths_selectivity.npy",
        allow_pickle=True,
    ).item()
    mean_1_ectopic = [np.mean(data["1_ectopic"][i]) for i in data["1_ectopic"].keys()]
    std_1_ectopic = [np.std(data["1_ectopic"][i]) for i in data["1_ectopic"].keys()]
    mean_2_ectopic = [np.mean(data["2_ectopic"][i]) for i in data["2_ectopic"].keys()]
    std_2_ectopic = [np.std(data["2_ectopic"][i]) for i in data["2_ectopic"].keys()]

    ax.plot(
        data["1_ectopic"].keys(),
        data["reference"],
        color="k",
        label="reference",
    )
    ax.plot(
        data["1_ectopic"].keys(),
        mean_1_ectopic,
        "-",
        label="1_ectopic",
        color="C3",
    )
    ax.errorbar(
        data["1_ectopic"].keys(),
        mean_1_ectopic,
        yerr=std_1_ectopic,
        fmt="o",
        color="C3",
        alpha=0.6,
    )
    ax.plot(
        data["2_ectopic"].keys(),
        mean_2_ectopic,
        "-",
        label="2_ectopic",
        color="C4",
    )
    ax.errorbar(
        data["2_ectopic"].keys(),
        mean_2_ectopic,
        yerr=std_2_ectopic,
        fmt="o",
        color="C4",
        alpha=0.6,
    )

    print("1 ectopic p-values")
    for i, j in enumerate(data["1_ectopic"].keys()):
        res = wilcoxon(
            np.array(data["1_ectopic"][j]) - data["reference"][i], alternative="less"
        )
        print(j, round(res.pvalue, 3))
        ax.text(
            x=j,
            y=1,
            s=convert_pvalue_to_stars(res.pvalue),
            va="center",
            ha="center",
            fontsize=10,
            color="C3",
        )
    print()
    print("2 ectopic p-values")
    for i, j in enumerate(data["2_ectopic"].keys()):
        res = wilcoxon(
            np.array(data["2_ectopic"][j]) - data["reference"][i], alternative="less"
        )
        print(j, round(res.pvalue, 3))
        ax.text(
            x=j,
            y=0.95,
            s=convert_pvalue_to_stars(res.pvalue),
            va="center",
            ha="center",
            fontsize=10,
            color="C4",
        )
    ax.legend(frameon=False)
    ax.set_xlabel("Sequence Length")
    ax.set_ylabel("Selectivity")
    ax.set_ylim(ymax=1.0)
    sns.despine()


def main():
    figure = plt.figure(figsize=(9, 8), constrained_layout=True)
    spec = figure.add_gridspec(6, 9)
    axa = figure.add_subplot(spec[0:2, 0:3])
    axb = figure.add_subplot(spec[0:2, 3:6])
    axc = figure.add_subplot(spec[0:2, 6:9])
    axd = figure.add_subplot(spec[2:4, 0:9])
    axe = figure.add_subplot(spec[4:6, 0:4])
    axf = figure.add_subplot(spec[4:6, 5:9])

    plot_panel_a(axa)
    plot_panel_b(axb)
    plot_panel_c(axc)
    plot_panel_d(axd)
    plot_panel_e(axe)
    plot_panel_f(axf)

    plt.figtext(0.02, 0.98, "A", fontsize=16, weight="bold")
    plt.figtext(0.35, 0.98, "B", fontsize=16, weight="bold")
    plt.figtext(0.66, 0.98, "C", fontsize=16, weight="bold")
    plt.figtext(0.05, 0.55, "D", fontsize=16, weight="bold")
    plt.figtext(0.02, 0.35, "E", fontsize=16, weight="bold")
    plt.figtext(0.55, 0.35, "F", fontsize=16, weight="bold")

    plt.savefig("Figure-5.png", bbox_inches="tight")
    plt.savefig("Figure-5.svg", bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
