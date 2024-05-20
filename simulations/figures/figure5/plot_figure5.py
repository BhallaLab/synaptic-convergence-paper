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


sns.set(style="ticks")
sns.set_context("paper")

DATA_PATH = "../../codes/sequences/chemical/data/"


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
    ax.legend(frameon=False)
    ax.set_xlabel("Position (\u03bcm)")
    ax.set_ylabel("A (a.u.)")
    sns.despine()


def plot_panel_c(ax):
    ax.imshow(mpimg.imread("ectopic_inputs.png"))
    ax.axis("off")


def plot_panel_d(ax):
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


def plot_panel_e(ax):
    """Plot selectivity in the absence of ectopic inputs, and in the presence of 1 or 2 ectopic inputs

    Parameters
    ----------
    ax : axis object
    """
    data = np.load(
        DATA_PATH + "seq_lengths_selectivity.npy",
        allow_pickle=True,
    ).item()
    ax.plot(data["1_ectopic"].keys(), data["reference"], color="k", label="reference")
    ax.plot(
        data["1_ectopic"].keys(),
        [np.mean(data["1_ectopic"][i]) for i in data["1_ectopic"].keys()],
        "-",
        label="1_ectopic",
    )
    ax.plot(
        data["2_ectopic"].keys(),
        [np.mean(data["2_ectopic"][i]) for i in data["2_ectopic"].keys()],
        "-",
        label="2_ectopic",
    )
    ax.legend(frameon=False)
    ax.set_xlabel("Sequence Length")
    ax.set_ylabel("Selectivity")
    sns.despine()


def main():
    figure = plt.figure(figsize=(7, 7), constrained_layout=True)
    spec = figure.add_gridspec(5, 8)
    axa = figure.add_subplot(spec[0:2, 0:3])
    axb = figure.add_subplot(spec[0:2, 3:8])
    axc = figure.add_subplot(spec[2:3, 0:8])
    axd = figure.add_subplot(spec[3:5, 0:4])
    axe = figure.add_subplot(spec[3:5, 4:8])
    plot_panel_a(axa)
    plot_panel_b(axb)
    plot_panel_c(axc)
    plot_panel_d(axd)
    plot_panel_e(axe)

    plt.figtext(0.02, 0.98, "A", fontsize=16, weight="bold")
    plt.figtext(0.36, 0.98, "B", fontsize=16, weight="bold")
    plt.figtext(0.15, 0.55, "C", fontsize=16, weight="bold")
    plt.figtext(0.02, 0.38, "D", fontsize=16, weight="bold")
    plt.figtext(0.52, 0.38, "E", fontsize=16, weight="bold")

    plt.savefig("Figure-5.png", bbox_inches="tight")
    plt.savefig("Figure-5.svg", bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
