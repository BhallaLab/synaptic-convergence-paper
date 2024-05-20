#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : plot_figure1.py
# Author            : Bhanu Priya S <bhanups@ncbs.res.in>
# Date              : 01.02.2023
# Last Modified Date: 01.02.2023
# Last Modified By  : Bhanu Priya S <bhanups@ncbs.res.in>

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.image as mpimg
import seaborn as sns
import cv2

sns.set(style="ticks")
sns.set_context("paper")


def plot_panel_a(ax):
    """Display stimuli schematic

    Parameters
    ----------
    ax : axis object
    """
    ax.imshow(mpimg.imread("Figure-1a.png"))
    ax.axis("off")


def plot_panel_b(ax):
    """Display convergence schematic

    Parameters
    ----------
    ax : axis object
    """
    ax.imshow(mpimg.imread("Figure-1b.png"))
    ax.axis("off")


def plot_panel_c(ax):
    """Display network architecture

    Parameters
    ----------
    ax : axis object
    """
    ax.imshow(mpimg.imread("Figure-1c.png"))
    ax.axis("off")


def main():
    figure = plt.figure(figsize=(7, 6), constrained_layout=True)
    spec = figure.add_gridspec(2, 2)
    axa = figure.add_subplot(spec[0, 0])
    axb = figure.add_subplot(spec[0, 1])
    axc = figure.add_subplot(spec[1, 0:2])

    plot_panel_a(axa)
    plot_panel_b(axb)
    plot_panel_c(axc)

    plt.figtext(0.04, 0.95, "A", fontsize=12, weight="bold")
    plt.figtext(0.55, 0.95, "B", fontsize=12, weight="bold")
    plt.figtext(0.04, 0.45, "C", fontsize=12, weight="bold")

    # Adding annotations and labels

    ax2 = plt.axes([0, 0, 1, 1], facecolor=(1, 1, 1, 0))

    x1, y1 = np.array([[0.750, 0.670], [0.267, 0.549]])
    line1 = mlines.Line2D(x1, y1, lw=1.0, color="k", alpha=1.0)
    x2, y2 = np.array([[0.759, 0.839], [0.267, 0.549]])
    line2 = mlines.Line2D(x2, y2, lw=1.0, color="k", alpha=1.0)
    ax2.add_line(line1)
    ax2.add_line(line2)

    ax2.text(0.27, 0.52, "Time", horizontalalignment="center")
    ax2.text(
        0.04,
        0.72,
        "Ensembles",
        horizontalalignment="center",
        rotation=90,
    )
    ax2.text(0.15, 0.95, "Grouped", horizontalalignment="center")
    ax2.text(0.35, 0.95, "Sequential", horizontalalignment="center")
    ax2.text(0.63, 0.58, "Group", horizontalalignment="center")
    ax2.text(0.86, 0.62, "Sequence", horizontalalignment="center")
    ax2.text(
        0.23,
        0.44,
        "Presynaptic population",
        horizontalalignment="center",
    )
    ax2.text(
        0.77,
        0.44,
        "Postsynaptic population",
        horizontalalignment="center",
    )
    ax2.text(
        0.5,
        0.28,
        "Random\nconnectivity\nprobability 'p'",
        horizontalalignment="center",
        verticalalignment="center",
    )
    ax2.text(0.48, 0.1, "Ensemble\nsize 'N'", horizontalalignment="center")
    ax2.axis("off")
    plt.savefig("Figure-1.svg", bbox_inches="tight")
    plt.savefig("Figure-1.png", bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
