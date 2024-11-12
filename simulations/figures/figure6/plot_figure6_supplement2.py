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


def plot_epsp_traces(
    file1, label1, file1_seq_indices, file2, label2, file2_seq_indices, ax
):
    """Plot soma response to inward, scrambled and outward sequences

    Parameters
    ----------
    file1 : str
        Name of the first file
    label1 : str
        Label for the first file
    file1_seq_indices : list
        List of sequence indices to plot from the first file for inward and outward sequences
    file2 : str
        Name of the second file
    label2 : str
        Label for the second file
    file2_seq_indices : list
        List of sequence indices to plot from the second file for inward and outward sequences
    ax : list
        List of axis objects
    """
    num_seq_dt = 9
    dt = 0.025

    first_file_responses = pd.read_csv(
        f"{DATA_PATH}{file1}",
        header=None,
        sep=" ",
    )
    for seq_dt_index in range(0, num_seq_dt):
        in_seq_col = file1_seq_indices[0] * num_seq_dt + seq_dt_index
        out_seq_col = file1_seq_indices[1] * num_seq_dt + seq_dt_index
        ax[0].plot(
            np.arange(first_file_responses.shape[0]) * dt,
            first_file_responses.iloc[:][in_seq_col],
            label=f"{seq_dt_index} ms",
        )
        ax[1].plot(
            np.arange(first_file_responses.shape[0]) * dt,
            first_file_responses.iloc[:][out_seq_col],
            label=f"{seq_dt_index} ms",
        )
    ax[0].set_ylim(-76, -62)
    ax[1].set_ylim(-76, -62)
    ax[0].set_ylabel("Soma potential (mV)")
    ax[0].text(
        0.4,
        0.9,
        s=f"Inward sequence - {label1}",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax[0].transAxes,
    )
    ax[1].text(
        0.5,
        0.9,
        s=f"Outward sequence - {label1}",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax[1].transAxes,
    )

    second_file_responses = pd.read_csv(
        f"{DATA_PATH}/{file2}",
        header=None,
        sep=" ",
    )
    for seq_dt_index in range(0, num_seq_dt):
        in_seq_col = file2_seq_indices[0] * num_seq_dt + seq_dt_index
        out_seq_col = file2_seq_indices[1] * num_seq_dt + seq_dt_index
        ax[2].plot(
            np.arange(second_file_responses.shape[0]) * dt,
            second_file_responses.iloc[:][in_seq_col],
            label=f"{seq_dt_index} ms",
        )
        ax[3].plot(
            np.arange(second_file_responses.shape[0]) * dt,
            second_file_responses.iloc[:][out_seq_col],
            label=f"{seq_dt_index} ms",
        )
    ax[2].set_ylim(-76, -62)
    ax[3].set_ylim(-76, -62)
    ax[2].text(
        0.5,
        0.9,
        s=f"Inward sequence - {label2}",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax[2].transAxes,
    )
    ax[3].text(
        0.5,
        0.9,
        s=f"Outward sequence - {label2}",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax[3].transAxes,
    )
    sns.despine()


def plot_dendritic_epsp_traces(files, labels, file_seq_indices, seq_dt_index, ax):
    """Plot soma response to inward, scrambled and outward sequences

    Parameters
    ----------
    file : list of str
        List of names of the data files
    label : list of str
        List of labels of the data files
    file_seq_indices : list of lists
        List of sequence indices for each file to plot for inward and outward sequences
    seq_dt_index : int
        Index of the seqence Dt to plot
    ax : list
        List of axis objects
    """
    num_seq_dt = 9
    dt = 0.025
    count = 0

    for filename, label in zip(files, labels):
        responses = pd.read_csv(
            f"{DATA_PATH}{filename}",
            header=None,
            sep=" ",
        )
        in_seq_col = file_seq_indices[count][0] * num_seq_dt + seq_dt_index
        out_seq_col = file_seq_indices[count][1] * num_seq_dt + seq_dt_index
        if "without_nmda" in filename:
            ax[2].plot(
                np.arange(responses.shape[0]) * dt,
                responses.iloc[:][in_seq_col],
                label=label,
            )
            ax[3].plot(
                np.arange(responses.shape[0]) * dt,
                responses.iloc[:][out_seq_col],
                # color="g",
                label=label,
            )
        else:
            ax[0].plot(
                np.arange(responses.shape[0]) * dt,
                responses.iloc[:][in_seq_col],
                label=label,
            )
            ax[1].plot(
                np.arange(responses.shape[0]) * dt,
                responses.iloc[:][out_seq_col],
                # color="g",
                label=label,
            )
        count += 1
    ax[0].set_ylim(-76, -25)
    ax[1].set_ylim(-76, -25)
    ax[2].set_ylim(-76, -25)
    ax[3].set_ylim(-76, -25)
    ax[0].set_ylabel("Dendritic potential (mV)")
    ax[0].text(
        0.5,
        0.9,
        s=f"Inward sequence - len 5\nAMPA + NMDA",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax[0].transAxes,
    )
    ax[1].text(
        0.5,
        0.9,
        s=f"Outward sequence - len 5\nAMPA + NMDA",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax[1].transAxes,
    )
    ax[2].text(
        0.5,
        0.9,
        s=f"Inward sequence - len 5\nAMPA",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax[2].transAxes,
    )
    ax[3].text(
        0.5,
        0.9,
        s=f"Outward sequence - len 5\nAMPA",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax[3].transAxes,
    )
    sns.despine()


def plot_input_velocity_vs_selectivity(
    file1, label1, file1_seq_indices, file2, label2, file2_seq_indices, ax
):
    """Plot selectivity for inward sequences as a function on input velocity

    Parameters
    ----------
    file1 : str
        Name of the first file
    label1 : str
        Label for the first file
    file1_seq_indices : list
        List of sequence indices to plot from the first file for inward and outward sequences
    file2 : str
        Name of the second file
    label2 : str
        Label for the second file
    file2_seq_indices : list
        List of sequence indices to plot from the second file for inward and outward sequences
    ax : axis object
    """
    L = 99.02251
    num_seq_dt = 9
    seq_dx = L / 9
    seq_dts = np.arange(num_seq_dt)
    first_file_responses = pd.read_csv(
        f"{DATA_PATH}{file1}",
        header=None,
        sep=" ",
    )
    selectivity = []
    for seq_dt_index, _ in enumerate(seq_dts):
        selectivity.append(
            first_file_responses.iloc[seq_dt_index][file1_seq_indices[0]]
            - first_file_responses.iloc[seq_dt_index][file1_seq_indices[1]]
        )
    ax.plot(
        [seq_dx / seq_dt for seq_dt in seq_dts],
        selectivity,
        marker="o",
        label=label1,
    )

    second_file_responses = pd.read_csv(
        f"{DATA_PATH}/{file2}",
        header=None,
        sep=" ",
    )

    selectivity = []
    for seq_dt_index, _ in enumerate(seq_dts):
        selectivity.append(
            second_file_responses.iloc[seq_dt_index][file2_seq_indices[0]]
            - second_file_responses.iloc[seq_dt_index][file2_seq_indices[1]]
        )
    ax.plot(
        [seq_dx / seq_dt for seq_dt in seq_dts],
        selectivity,
        marker="o",
        label=label2,
    )
    ax.set_xlabel("Input velocity (\u03bcm/ms)")
    ax.set_ylim(0, 3)
    sns.despine()


def main():
    figure = plt.figure(figsize=(15, 14), constrained_layout=True)
    spec = figure.add_gridspec(5, 12)
    axa = figure.add_subplot(spec[0:1, 0:3])
    axb = figure.add_subplot(spec[0:1, 3:6])
    axc = figure.add_subplot(spec[0:1, 6:9])
    axd = figure.add_subplot(spec[0:1, 9:12])
    axe = figure.add_subplot(spec[1:2, 0:3])
    axf = figure.add_subplot(spec[1:2, 3:6])
    axg = figure.add_subplot(spec[1:2, 6:9])
    axh = figure.add_subplot(spec[1:2, 9:12])
    axi = figure.add_subplot(spec[2:3, 0:3])
    axj = figure.add_subplot(spec[2:3, 3:6])
    axk = figure.add_subplot(spec[2:3, 6:9])
    axl = figure.add_subplot(spec[2:3, 9:12])
    axm = figure.add_subplot(spec[3:4, 0:4])
    axn = figure.add_subplot(spec[3:4, 4:8])
    axo = figure.add_subplot(spec[3:4, 8:12])
    axp = figure.add_subplot(spec[4:5, 0:3])
    axq = figure.add_subplot(spec[4:5, 3:6])
    axr = figure.add_subplot(spec[4:5, 6:9])
    axs = figure.add_subplot(spec[4:5, 9:12])

    plot_epsp_traces(
        "tempV_5_panelC.dat",
        "len 5\nAMPA + NMDA",
        [0, 119],
        "tempV_5_panelC_without_nmda.dat",
        "len 5\nAMPA",
        [0, 119],
        [axa, axb, axc, axd],
    )
    axa.legend(frameon=False)
    plot_input_velocity_vs_selectivity(
        "tempM_5_panelC.dat",
        "AMPA + NMDA",
        [0, 119],
        "tempM_5_panelC_without_nmda.dat",
        "AMPA",
        [0, 119],
        axm,
    )
    axm.set_ylabel("$V_{IN}$ - $V_{OUT}$ (mV)")
    axm.set_title("Seq len 5")
    axm.legend(frameon=False)
    plot_epsp_traces(
        "tempV_9_panelC.dat",
        "len 9\nAMPA + NMDA",
        [0, 120],
        "tempV_9_panelC_without_nmda.dat",
        "len 9\nAMPA",
        [0, 120],
        [axe, axf, axg, axh],
    )
    plot_input_velocity_vs_selectivity(
        "tempM_9_panelC.dat",
        "AMPA + NMDA",
        [0, 120],
        "tempM_9_panelC_without_nmda.dat",
        "AMPA",
        [0, 120],
        axn,
    )
    axn.set_title("Seq len 9")
    plot_epsp_traces(
        "branco-model/tempV.dat",
        "len 9\nAMPA + NMDA, Branco et al",
        [0, 1],
        "branco-model/tempV_without_nmda.dat",
        "len 9\nAMPA, Branco et al",
        [0, 1],
        [axi, axj, axk, axl],
    )
    plot_input_velocity_vs_selectivity(
        "branco-model/tempM.dat",
        "AMPA + NMDA",
        [0, 1],
        "branco-model/tempM_without_nmda.dat",
        "AMPA",
        [0, 1],
        axo,
    )
    axo.set_title("Seq len 9, Branco et al.")
    plot_dendritic_epsp_traces(
        [
            "tempdendV_5_panelC_25.dat",
            "tempdendV_5_panelC_50.dat",
            "tempdendV_5_panelC_75.dat",
            "tempdendV_5_panelC_without_nmda_25.dat",
            "tempdendV_5_panelC_without_nmda_50.dat",
            "tempdendV_5_panelC_without_nmda_75.dat",
        ],
        ["0.25", "0.5", "0.75", "0.25", "0.5", "0.75"],
        [[0, 119], [0, 119], [0, 119], [0, 119], [0, 119], [0, 119]],
        3,
        [axp, axq, axr, axs],
    )
    axp.legend(loc="right", frameon=False)
    for cx in [axi, axj, axk, axl, axp, axq, axr, axs]:
        cx.set_xlabel("Time (ms)")

    plt.figtext(0.00, 0.98, "A", fontsize=16, weight="bold")
    plt.figtext(0.24, 0.98, "B", fontsize=16, weight="bold")
    plt.figtext(0.49, 0.98, "C", fontsize=16, weight="bold")
    plt.figtext(0.76, 0.98, "D", fontsize=16, weight="bold")
    plt.figtext(0.00, 0.79, "E", fontsize=16, weight="bold")
    plt.figtext(0.24, 0.79, "F", fontsize=16, weight="bold")
    plt.figtext(0.49, 0.79, "G", fontsize=16, weight="bold")
    plt.figtext(0.76, 0.79, "H", fontsize=16, weight="bold")
    plt.figtext(0.00, 0.59, "I", fontsize=16, weight="bold")
    plt.figtext(0.24, 0.59, "J", fontsize=16, weight="bold")
    plt.figtext(0.49, 0.59, "K", fontsize=16, weight="bold")
    plt.figtext(0.76, 0.59, "L", fontsize=16, weight="bold")
    plt.figtext(0.00, 0.38, "M", fontsize=16, weight="bold")
    plt.figtext(0.33, 0.38, "N", fontsize=16, weight="bold")
    plt.figtext(0.68, 0.38, "O", fontsize=16, weight="bold")
    plt.figtext(0.00, 0.17, "P", fontsize=16, weight="bold")
    plt.figtext(0.24, 0.17, "Q", fontsize=16, weight="bold")
    plt.figtext(0.49, 0.17, "R", fontsize=16, weight="bold")
    plt.figtext(0.76, 0.17, "S", fontsize=16, weight="bold")

    plt.savefig("Figure-6-supplement2.png", bbox_inches="tight")
    plt.savefig("Figure-6-supplement2.svg", bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
