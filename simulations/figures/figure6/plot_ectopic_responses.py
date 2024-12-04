#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : plot_ectopic_responses.py
# Author            : Bhanu Priya S <bhanups@ncbs.res.in>
# Date              : 01.02.2023
# Last Modified Date: 21.09.2023
# Last Modified By  : Bhanu Priya S <bhanups@ncbs.res.in>

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

sns.set(style="ticks")
sns.set_context("paper")

DATA_PATH = "../../codes/sequences/electrical/"


def plotPanelA(ax):
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
        -0.1,
        "proximal",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    ax.text(
        0.9,
        -0.1,
        "distal",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    sns.despine()


def plotPanelB(ax, bx):
    """Plot vmax in the absence and presence of ectopic inputs

    Parameters
    ----------
    ax : axis object
    bx : axis object
    """
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
    colors = ["b", "orange", "g", "r"]

    ax = sns.violinplot(
        responses,
        x="ect_pos",
        y="vmax",
        hue="ect_time",
        palette=sns.color_palette(),
        inner=None,
        ax=ax,
    )
    bx = sns.swarmplot(
        responses,
        x="ect_pos",
        y="vmax",
        dodge=True,
        hue="ect_time",
        palette=sns.color_palette(),
        ax=bx,
        size=1,
    )
    plt.setp(ax.collections, alpha=0.7)
    plt.setp(ax.collections, edgecolor=(1, 1, 1))
    i = 0
    for ect_info, group in responses.groupby(["ect_pos", "ect_time"]):
        ect_pos, ect_time = ect_info
        ax.scatter(
            int(i / 4) + (-0.3 + (i % 4) * 0.2),
            group.loc[group["seq_id"] == 0, "vmax"],
            color=colors[i % 4],
            marker="v",
        )
        bx.scatter(
            int(i / 4) + (-0.3 + (i % 4) * 0.2),
            group.loc[group["seq_id"] == 0, "vmax"],
            color=colors[i % 4],
            marker="v",
        )
        i += 1

    ax.legend(
        handles=ax.legend_.legendHandles,
        labels=["t_start", "t_mid", "t_end", "reference"],
        frameon=False,
    )
    ax.set_ylabel("EPSP (mV)")
    ax.set_xlabel("")
    bx.legend(
        handles=bx.legend_.legendHandles,
        labels=["t_start", "t_mid", "t_end", "reference"],
        frameon=False,
    )
    bx.set_xlabel("Relative position of ectopic input on dendrite")
    bx.set_ylabel("EPSP (mV)")
    sns.despine()


def main():
    figure = plt.figure(figsize=(15, 9), constrained_layout=True)
    spec = figure.add_gridspec(3, 1)
    axa = figure.add_subplot(spec[0, 0])
    axb = figure.add_subplot(spec[1, 0])
    axc = figure.add_subplot(spec[2, 0])
    plotPanelA(axa)
    plotPanelB(axb, axc)

    plt.figtext(0.01, 0.98, "A", fontsize=16, weight="bold")
    plt.figtext(0.01, 0.68, "B", fontsize=16, weight="bold")
    plt.figtext(0.01, 0.33, "C", fontsize=16, weight="bold")

    plt.savefig("ectopic_responses.svg", bbox_inches="tight")
    plt.savefig("ectopic_responses.png", bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
