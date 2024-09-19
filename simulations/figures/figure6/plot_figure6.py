#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : plot_figure7.py
# Author            : Bhanu Priya S <bhanups@ncbs.res.in>
# Date              : 01.02.2023
# Last Modified Date: 10.01.2024
# Last Modified By  : Bhanu Priya S <bhanups@ncbs.res.in>

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import seaborn as sns
import pandas as pd
from scipy.stats import wilcoxon

sys.path.append("../../codes/sequences/")
import sequence_params as prm

sns.set(style="ticks")
sns.set_context("paper")

DATA_PATH = "../../codes/sequences/electrical/"

MIN_M = prm.MIN_M
MAX_M = prm.MAX_M


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


def plotPanelA(ax):
    """Display schematic of bistable-switch model

    Parameters
    ----------
    ax : axis object
    """
    ax.imshow(mpimg.imread("Figure-6a.png"))
    ax.axis("off")


def plotPanelB(ax):
    """Plot soma response to inward, scrambled and outward sequences

    Parameters
    ----------
    ax : axis object
    """
    responses = pd.read_csv(
        DATA_PATH + "tempV_5_panelC.dat",
        header=None,
        sep=" ",
    )
    num_seq_dt = 9
    seq_dt_index = 3
    dt = 0.025
    in_seq_col = 0 * num_seq_dt + seq_dt_index
    scram_seq_col = 58 * num_seq_dt + seq_dt_index
    out_seq_col = 119 * num_seq_dt + seq_dt_index
    sequences = {
        "inward": in_seq_col,
        "scrambled": scram_seq_col,
        "outward": out_seq_col,
    }

    for label, seq_col in sequences.items():
        ax.plot(
            np.arange(responses.shape[0]) * dt,
            responses.iloc[:][seq_col],
            label=label,
        )
    ax.legend(frameon=False)
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Soma potential (mV)")
    sns.despine()


def plotPanelC(ax):
    """Plot selectivity for inward sequences as a function on input velocity

    Parameters
    ----------
    ax : axis object
    """
    L = 99.02251
    num_seq_dt = 9
    seq_dx = L / 9
    seq_dts = np.arange(num_seq_dt)
    for M in [3, 5, 8]:
        responses = pd.read_csv(
            f"{DATA_PATH}tempM_{M}_panelC.dat",
            header=None,
            sep=" ",
        )
        selectivity = []
        for seq_dt_index, _ in enumerate(seq_dts):
            selectivity.append(
                (
                    responses.iloc[seq_dt_index][0]
                    - np.mean(responses.iloc[seq_dt_index][:])
                )
                / np.max(responses.iloc[seq_dt_index][:])
            )
        ax.plot(
            [seq_dx / seq_dt for seq_dt in seq_dts],
            selectivity,
            marker="o",
            label=M,
        )
    ax.legend(frameon=False)
    ax.set_xlabel("Input velocity (\u03bcm/ms)")
    ax.set_ylabel("Sequence selectivity (a.u.)")
    sns.despine()


def plotPanelD(ax):
    """Plot selectivity in the absence and presence of ectopic inputs

    Parameters
    ----------
    ax : axis object
    """
    L = 99.02251
    M = 5
    responses = pd.read_csv(
        DATA_PATH + "tempM_panelD.dat",
        header=None,
        sep=" ",
        usecols=range(4),
        dtype={
            0: "int64",
            1: "float64",
            2: "int64",
            3: "float64",
        },
    )
    responses.columns = ["seq_id", "ect_pos", "ect_time", "vmax"]
    ect_time_labels = {
        0: "t_start",
        2: "t_mid",
        4: "t_end",
        10000: "reference",
    }
    synapse_locs = pd.read_csv(
        DATA_PATH + f"synapse_loc_{M}.dat",
        header=None,
        sep=" ",
        skiprows=[0],
        usecols=[1],
        dtype={
            0: "float64",
        },
    )
    for ect_time, time_group in responses.groupby("ect_time"):
        selectivity = []
        ectopic_positions = []
        for ect_pos, pos_group in time_group.groupby("ect_pos"):
            sorted_group = pos_group.sort_values(by=["seq_id"])
            sel = (
                sorted_group.iloc[0]["vmax"] - np.mean(sorted_group.iloc[:]["vmax"])
            ) / np.max(sorted_group.iloc[:]["vmax"])
            ectopic_positions.append(ect_pos * L)
            selectivity.append(sel)
        if ect_time == 10000:
            ax.plot(
                ectopic_positions,
                selectivity,
                marker="o",
                markersize=3,
                label=ect_time_labels[ect_time],
                color="k",
            )
        else:
            ax.plot(
                ectopic_positions,
                selectivity,
                marker="o",
                markersize=3,
                label=ect_time_labels[ect_time],
            )
    ax.axvspan(
        synapse_locs.iat[0, 0] * L,
        synapse_locs.iat[-1, 0] * L,
        color="c",
        alpha=0.5,
        lw=0,
    )
    ax.legend(frameon=False)
    ax.set_xlabel("Position of ectopic input (\u03bcm)")
    ax.set_ylabel("Sequence selectivity (a.u.)")
    ax.set_ylim(0)
    ax.text(
        0.1,
        1.1,
        "proximal",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.text(
        0.9,
        1.1,
        "distal",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    sns.despine()


def plotPanelE(ax):
    """Plot selectivity in the absence of ectopic inputs, and in the presence of 1 or 2 ectopic inputs

    Parameters
    ----------
    ax : axis object
    """
    selectivity = {}
    selectivity["reference"] = []
    selectivity["1_ectopic"] = []
    std_selectivity_1_ectopic = []
    for i, M in enumerate(range(MIN_M, MAX_M + 1)):
        responses = pd.read_csv(
            DATA_PATH + f"tempM_{M}_panelE.dat",
            header=None,
            sep=" ",
            usecols=range(4),
            dtype={
                0: "int64",
                1: "float64",
                2: "float64",
                3: "float64",
            },
        )
        responses.columns = ["seq_id", "ect_pos", "ect_time", "vmax"]
        ect_times = pd.read_csv(
            DATA_PATH + f"ectopic_input_time_{M}.dat",
            header=None,
            sep=" ",
            skiprows=[0],
            dtype={
                0: "float64",
            },
        )
        # Based on the order in which they have been defined
        for ect_type in ["reference", "1_ectopic"]:
            if ect_type == "reference":
                data = responses[responses["ect_time"] > 1000]
                seq_response = np.mean(data.loc[data["seq_id"] == 0, "vmax"])
                mean_response = np.mean(data["vmax"])
                max_response = np.max(data["vmax"])
                selectivity[ect_type].append(
                    (seq_response - mean_response) / max_response
                )
            else:
                data = responses[responses["ect_time"] < 1000]
                sel_group = []
                for ect_time, time_group in data.groupby("ect_time"):
                    for ect_pos, pos_group in time_group.groupby("ect_pos"):
                        sorted_group = pos_group.sort_values(by=["seq_id"])
                        sel = (
                            sorted_group.iloc[0]["vmax"]
                            - np.mean(sorted_group.iloc[:]["vmax"])
                        ) / np.max(sorted_group.iloc[:]["vmax"])
                        sel_group.append(sel)
                selectivity[ect_type].append(np.mean(sel_group))
                std_selectivity_1_ectopic.append(np.std(sel_group))
                res = wilcoxon(
                    np.array(sel_group) - selectivity["reference"][i],
                    alternative="less",
                )
                print(M, round(res.pvalue, 3))
                ax.text(
                    x=M,
                    y=0.1,
                    s=convert_pvalue_to_stars(res.pvalue),
                    va="center",
                    ha="center",
                    fontsize=10,
                    color="C3",
                )

    ax.plot(
        range(MIN_M, MAX_M + 1), selectivity["reference"], color="k", label="reference"
    )
    ax.plot(
        range(MIN_M, MAX_M + 1), selectivity["1_ectopic"], color="C3", label="1_ectopic"
    )
    ax.errorbar(
        range(MIN_M, MAX_M + 1),
        selectivity["1_ectopic"],
        yerr=std_selectivity_1_ectopic,
        color="r",
        fmt="o",
        alpha=0.5,
    )
    ax.legend(frameon=False)
    ax.set_xlabel("Sequence length")
    ax.set_ylabel("Sequence selectivity (a.u.)")
    ax.set_ylim(0)
    sns.despine()


def main():
    figure = plt.figure(figsize=(7, 7), constrained_layout=True)
    spec = figure.add_gridspec(6, 8)
    axa = figure.add_subplot(spec[0:4, 0:4])
    axb = figure.add_subplot(spec[0:2, 4:8])
    axc = figure.add_subplot(spec[2:4, 4:8])
    axd = figure.add_subplot(spec[4:6, 0:4])
    axe = figure.add_subplot(spec[4:6, 4:8])
    plotPanelA(axa)
    plotPanelB(axb)
    plotPanelC(axc)
    plotPanelD(axd)
    plotPanelE(axe)

    plt.figtext(0.01, 0.98, "A", fontsize=16, weight="bold")
    plt.figtext(0.51, 0.98, "B", fontsize=16, weight="bold")
    plt.figtext(0.51, 0.67, "C", fontsize=16, weight="bold")
    plt.figtext(0.01, 0.33, "D", fontsize=16, weight="bold")
    plt.figtext(0.51, 0.33, "E", fontsize=16, weight="bold")

    plt.savefig("Figure-6.png", bbox_inches="tight")
    plt.savefig("Figure-6.svg", bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
