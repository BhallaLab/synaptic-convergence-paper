#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : plot_figure2_param_screen.py
# Author            : Bhanu Priya S <bhanups@ncbs.res.in>
# Date              : 09.11.2023
# Last Modified Date: 09.11.2023
# Last Modified By  : Bhanu Priya S <bhanups@ncbs.res.in>

"""Plot Figure-2
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib import cm
from matplotlib.colors import LogNorm
import seaborn as sns
import pandas as pd
from scipy.linalg import expm
import matplotlib.ticker as mticker
from scipy.stats import poisson


sys.path.append("../../codes/groups/")
import group_params as prm

sns.set(style="ticks")
sns.set_context("paper")

DATA_PATH = "../../codes/groups/data/"
CAM_MODEL_DATA_PATH = "../../codes/groups/cam/data/"

networks = [
    prm.hippo_chem,
    prm.hippo_cicr,
    prm.hippo_elec,
    prm.cortex_chem,
    prm.cortex_cicr,
    prm.cortex_elec,
]
colors = sns.color_palette("muted", n_colors=len(networks))

MIN_M = prm.MIN_M
MAX_M = prm.MAX_M
T_POP_POST = 4000
NUM_RUNS = 100


def plot_conn_prob_vs_p_group(ax, net, m):
    """Plot connection probability vs probability of group

    Parameters
    ----------
    ax : axis object
    net : dict
        Dictionary with network configuration parameters
    m : int
        Group size
    """
    p = np.logspace(-2, 0, num=7)
    p_e = net["ensemble_participation_prob"]
    N = net["ensemble_size"]
    L = net["dendrite_length"]
    Z = net["zone_length"]
    t_pop_pre = net["population_size"]
    R = net["background_rate"]
    D = net["input_duration"]

    # True connected group
    p_total_exp_fit = []
    for p_i in p:
        sigma = L / int(p_i * t_pop_pre)
        q = N / t_pop_pre
        x = q * L / sigma
        y = L / Z
        A = np.zeros((m + 1, m + 1))
        for j in range(m):
            A[j, j] = -(m - j) * x - j * y
            A[j + 1, j] = (m - j) * x
        for j in range(m - 1):
            A[j, j + 1] = (j + 1) * y
        B = expm(A)
        p_total_exp_fit.append(B[m, 0])
    # ax.plot(p, p_total_exp_fit, "-o", c="b", label="true connected")
    # c = 0.7384796
    c = 2.56141683
    egrp_exp_Ei = np.array([p_i * N * Z / L for p_i in p])
    p_grp_exp = np.power(1 - np.exp(-egrp_exp_Ei), m)
    p_total_exp_1 = 1 - np.power((1 - p_grp_exp), c * int(L / Z))
    ax.plot(p, p_total_exp_1, "-o", c="b", label="true connected", alpha=0.7)

    # True active group
    p_total_exp_fit = []
    for p_i in p:
        sigma = L / int(p_i * t_pop_pre)
        q = p_e * N / t_pop_pre
        x = q * L / sigma
        y = L / Z
        A = np.zeros((m + 1, m + 1))
        for j in range(m):
            A[j, j] = -(m - j) * x - j * y
            A[j + 1, j] = (m - j) * x
        for j in range(m - 1):
            A[j, j + 1] = (j + 1) * y
        B = expm(A)
        p_total_exp_fit.append(B[m, 0])
    # ax.plot(p, p_total_exp_fit, ":o", c="b", label="true active")
    # c = 0.74330444
    c = 2.57658637
    egrp_exp_Ei = np.array([p_i * p_e * N * Z / L for p_i in p])
    p_grp_exp = np.power(1 - np.exp(-egrp_exp_Ei), m)
    p_total_exp_1 = 1 - np.power((1 - p_grp_exp), c * int(L / Z))
    ax.plot(p, p_total_exp_1, ":o", c="b", label="true active", alpha=0.7)

    # Degenerate connected group
    p_total_exp_fit = []
    for p_i in p:
        sigma = L / int(p_i * t_pop_pre)
        q = (m * N) / t_pop_pre
        x = q * L / sigma
        y = L / Z
        A = np.zeros((m + 1, m + 1))
        for j in range(m):
            A[j, j] = -x - j * y
            A[j + 1, j] = x
        for j in range(m - 1):
            A[j, j + 1] = (j + 1) * y
        B = expm(A)
        p_total_exp_fit.append(B[m, 0])
    # ax.plot(p, p_total_exp_fit, "-s", c="c", label="degenerate connected")
    # c = 0.66733867
    c = 2.62999896
    p_total_exp_1 = []
    for p_i in p:
        ealt_all_exp_EM = p_i * N * m * Z / L
        p_alt_exp = 1 - np.sum(
            [poisson.pmf(k=j, mu=ealt_all_exp_EM) for j in range(0, m)]
        )
        p_total_exp_1.append(1 - np.power((1 - p_alt_exp), c * int(L / Z)))
    ax.plot(p, p_total_exp_1, "-s", c="c", label="degenerate connected", alpha=0.7)

    # Degenerate active group
    p_total_exp_fit = []
    for p_i in p:
        sigma = L / int(p_i * t_pop_pre)
        q = (p_e * m * N) / t_pop_pre
        x = q * L / sigma
        y = L / Z
        A = np.zeros((m + 1, m + 1))
        for j in range(m):
            A[j, j] = -x - j * y
            A[j + 1, j] = x
        for j in range(m - 1):
            A[j, j + 1] = (j + 1) * y
        B = expm(A)
        p_total_exp_fit.append(B[m, 0])
    # ax.plot(p, p_total_exp_fit, ":s", c="c", label="degenerate active")
    # c = 0.67088243
    c = 2.63552124
    p_total_exp_1 = []
    for p_i in p:
        ealt_all_exp_EM = p_i * p_e * N * m * Z / L
        p_alt_exp = 1 - np.sum(
            [poisson.pmf(k=j, mu=ealt_all_exp_EM) for j in range(0, m)]
        )
        p_total_exp_1.append(1 - np.power((1 - p_alt_exp), c * int(L / Z)))
    ax.plot(p, p_total_exp_1, ":s", c="c", label="degenerate active", alpha=0.7)

    # Noise group
    p_total_exp_fit = []
    for p_i in p:
        sigma = L / int(p_i * t_pop_pre)
        q = ((t_pop_pre - m * N) * (1 - np.exp(-R * D))) / t_pop_pre
        x = q * L / sigma
        y = L / Z
        A = np.zeros((m + 1, m + 1))
        for j in range(m):
            A[j, j] = -x - j * y
            A[j + 1, j] = x
        for j in range(m - 1):
            A[j, j + 1] = (j + 1) * y
        B = expm(A)
        p_total_exp_fit.append(B[m, 0])
    # ax.plot(p, p_total_exp_fit, ":+", c="r", label="noise")
    # c = 0.39755312
    c = 2.30756123
    p_total_exp_1 = []
    for p_i in p:
        sigma = L / int(p_i * t_pop_pre)
        ENoise_exp = R * D * Z * ((1 / sigma) - (p_i * N * m / L))
        p_noise_exp = 1 - np.sum([poisson.pmf(k=j, mu=ENoise_exp) for j in range(0, m)])
        p_total_exp_1.append(1 - np.power((1 - p_noise_exp), c * int(L / Z)))
    ax.plot(p, p_total_exp_1, ":+", c="r", label="noise", alpha=0.7)

    # Any group
    p_total_exp_fit = []
    for p_i in p:
        sigma = L / int(p_i * t_pop_pre)
        q = (p_e * m * N + (t_pop_pre - m * N) * (1 - np.exp(-R * D))) / t_pop_pre
        x = q * L / sigma
        y = L / Z
        A = np.zeros((m + 1, m + 1))
        for j in range(m):
            A[j, j] = -x - j * y
            A[j + 1, j] = x
        for j in range(m - 1):
            A[j, j + 1] = (j + 1) * y
        B = expm(A)
        p_total_exp_fit.append(B[m, 0])
    # ax.plot(p, p_total_exp_fit, ":x", c="k", label="any")
    # c = 0.32892748
    c = 2.10499296
    p_total_exp_1 = []
    for p_i in p:
        sigma = L / int(p_i * t_pop_pre)
        EAll_exp = R * D * Z * ((1 / sigma) - (p_i * N * m / L)) + (
            p_i * p_e * N * m * Z / L
        )
        p_any_exp = 1 - np.sum([poisson.pmf(k=j, mu=EAll_exp) for j in range(0, m)])
        p_total_exp_1.append(1 - np.power((1 - p_any_exp), c * int(L / Z)))
    ax.plot(p, p_total_exp_1, ":x", c="k", label="any", alpha=0.7)

    sns.despine()
    ax.set_xscale("log", nonpositive="mask")
    ax.set_xlabel("Connection probability")
    ax.set_yscale("log", nonpositive="mask")
    ax.set_ylim([1e-6, 1e0])
    ax.set_ylabel("p_group")
    ax.legend(title="Legend for A,B,C,D", frameon=False)


def plot_ensemble_participation_prob_vs_p_group(ax, net, m):
    """Plot ensemble participation probability vs probability of group

    Parameters
    ----------
    ax : axis object
    net : dict
        Dictionary with network configuration parameters
    m : int
        Group size
    """
    p = net["connection_prob"]
    p_e = np.linspace(0, 1, 11)
    N = net["ensemble_size"]
    L = net["dendrite_length"]
    Z = net["zone_length"]
    t_pop_pre = net["population_size"]
    R = net["background_rate"]
    D = net["input_duration"]

    # True connected group
    # c = 0.7384796
    c = 2.56141683
    p_total_exp_1 = []
    for p_i in p_e:
        egrp_exp_Ei = p * N * Z / L
        p_grp_exp = np.power(1 - np.exp(-egrp_exp_Ei), m)
        p_total_exp_1.append(1 - np.power((1 - p_grp_exp), c * int(L / Z)))
    ax.plot(p_e, p_total_exp_1, "-o", c="b", label="true connected", alpha=0.7)

    # True active group
    # c = 0.74330444
    c = 2.57658637
    p_total_exp_1 = []
    for p_i in p_e:
        egrp_exp_Ei = p * p_i * N * Z / L
        p_grp_exp = np.power(1 - np.exp(-egrp_exp_Ei), m)
        p_total_exp_1.append(1 - np.power((1 - p_grp_exp), c * int(L / Z)))
    ax.plot(p_e, p_total_exp_1, ":o", c="b", label="true active", alpha=0.7)

    # Degenerate connected group
    # c = 0.66733867
    c = 2.62999896
    p_total_exp_1 = []
    for p_i in p_e:
        ealt_all_exp_EM = p * N * m * Z / L
        p_alt_exp = 1 - np.sum(
            [poisson.pmf(k=j, mu=ealt_all_exp_EM) for j in range(0, m)]
        )
        p_total_exp_1.append(1 - np.power((1 - p_alt_exp), c * int(L / Z)))
    ax.plot(p_e, p_total_exp_1, "-s", c="c", label="degenerate connected", alpha=0.7)

    # Degenerate active group
    # c = 0.67088243
    c = 2.63552124
    p_total_exp_1 = []
    for p_i in p_e:
        ealt_all_exp_EM = p * p_i * N * m * Z / L
        p_alt_exp = 1 - np.sum(
            [poisson.pmf(k=j, mu=ealt_all_exp_EM) for j in range(0, m)]
        )
        p_total_exp_1.append(1 - np.power((1 - p_alt_exp), c * int(L / Z)))
    ax.plot(p_e, p_total_exp_1, ":s", c="c", label="degenerate active", alpha=0.7)

    # Noise group
    # c = 0.39755312
    c = 2.30756123
    p_total_exp_1 = []
    for p_i in p_e:
        sigma = L / int(p * t_pop_pre)
        ENoise_exp = R * D * Z * ((1 / sigma) - (p * N * m / L))
        p_noise_exp = 1 - np.sum([poisson.pmf(k=j, mu=ENoise_exp) for j in range(0, m)])
        p_total_exp_1.append(1 - np.power((1 - p_noise_exp), c * int(L / Z)))
    ax.plot(p_e, p_total_exp_1, ":+", c="r", label="noise", alpha=0.7)

    # Any group
    # c = 0.32892748
    c = 2.10499296
    p_total_exp_1 = []
    for p_i in p_e:
        sigma = L / int(p * t_pop_pre)
        EAll_exp = R * D * Z * ((1 / sigma) - (p * N * m / L)) + (
            p * p_i * N * m * Z / L
        )
        p_any_exp = 1 - np.sum([poisson.pmf(k=j, mu=EAll_exp) for j in range(0, m)])
        p_total_exp_1.append(1 - np.power((1 - p_any_exp), c * int(L / Z)))
    ax.plot(p_e, p_total_exp_1, ":x", c="k", label="any", alpha=0.7)

    # # True connected group
    # p_total_exp_fit = []
    # for p_i in p_e:
    #     sigma = L / int(p * t_pop_pre)
    #     q = N / t_pop_pre
    #     x = q * L / sigma
    #     y = L / Z
    #     A = np.zeros((m + 1, m + 1))
    #     for j in range(m):
    #         A[j, j] = -(m - j) * x - j * y
    #         A[j + 1, j] = (m - j) * x
    #     for j in range(m - 1):
    #         A[j, j + 1] = (j + 1) * y
    #     B = expm(A)
    #     p_total_exp_fit.append(B[m, 0])
    # ax.plot(p_e, p_total_exp_fit, "-o", c="b", label="true connected")
    #
    # # True active group
    # p_total_exp_fit = []
    # for p_i in p_e:
    #     sigma = L / int(p * t_pop_pre)
    #     q = p_i * N / t_pop_pre
    #     x = q * L / sigma
    #     y = L / Z
    #     A = np.zeros((m + 1, m + 1))
    #     for j in range(m):
    #         A[j, j] = -(m - j) * x - j * y
    #         A[j + 1, j] = (m - j) * x
    #     for j in range(m - 1):
    #         A[j, j + 1] = (j + 1) * y
    #     B = expm(A)
    #     p_total_exp_fit.append(B[m, 0])
    # ax.plot(p_e, p_total_exp_fit, ":o", c="b", label="true active")
    #
    # # Degenerate connected group
    # p_total_exp_fit = []
    # for p_i in p_e:
    #     sigma = L / int(p * t_pop_pre)
    #     q = (m * N) / t_pop_pre
    #     x = q * L / sigma
    #     y = L / Z
    #     A = np.zeros((m + 1, m + 1))
    #     for j in range(m):
    #         A[j, j] = -x - j * y
    #         A[j + 1, j] = x
    #     for j in range(m - 1):
    #         A[j, j + 1] = (j + 1) * y
    #     B = expm(A)
    #     p_total_exp_fit.append(B[m, 0])
    # ax.plot(p_e, p_total_exp_fit, "-s", c="c", label="degenerate connected")
    #
    # # Degenerate active group
    # p_total_exp_fit = []
    # for p_i in p_e:
    #     sigma = L / int(p * t_pop_pre)
    #     q = (p_i * m * N) / t_pop_pre
    #     x = q * L / sigma
    #     y = L / Z
    #     A = np.zeros((m + 1, m + 1))
    #     for j in range(m):
    #         A[j, j] = -x - j * y
    #         A[j + 1, j] = x
    #     for j in range(m - 1):
    #         A[j, j + 1] = (j + 1) * y
    #     B = expm(A)
    #     p_total_exp_fit.append(B[m, 0])
    # ax.plot(p_e, p_total_exp_fit, ":s", c="c", label="degenerate active")
    #
    # # Noise group
    # p_total_exp_fit = []
    # for p_i in p_e:
    #     sigma = L / int(p * t_pop_pre)
    #     q = ((t_pop_pre - m * N) * (1 - np.exp(-R * D))) / t_pop_pre
    #     x = q * L / sigma
    #     y = L / Z
    #     A = np.zeros((m + 1, m + 1))
    #     for j in range(m):
    #         A[j, j] = -x - j * y
    #         A[j + 1, j] = x
    #     for j in range(m - 1):
    #         A[j, j + 1] = (j + 1) * y
    #     B = expm(A)
    #     p_total_exp_fit.append(B[m, 0])
    # ax.plot(p_e, p_total_exp_fit, ":+", c="r", label="noise")
    #
    # # Any group
    # p_total_exp_fit = []
    # for p_i in p_e:
    #     sigma = L / int(p * t_pop_pre)
    #     q = (p_i * m * N + (t_pop_pre - m * N) * (1 - np.exp(-R * D))) / t_pop_pre
    #     x = q * L / sigma
    #     y = L / Z
    #     A = np.zeros((m + 1, m + 1))
    #     for j in range(m):
    #         A[j, j] = -x - j * y
    #         A[j + 1, j] = x
    #     for j in range(m - 1):
    #         A[j, j + 1] = (j + 1) * y
    #     B = expm(A)
    #     p_total_exp_fit.append(B[m, 0])
    # ax.plot(p_e, p_total_exp_fit, ":x", c="k", label="any")

    sns.despine()
    # ax.set_xscale("log", nonpositive="mask")
    ax.set_xlabel("Ensemble participation\nprobability p_e")
    ax.set_yscale("log", nonpositive="mask")
    ax.set_ylim([1e-6, 1e0])
    ax.set_ylabel("p_group")
    # ax.legend(frameon=False)


def plot_zone_length_vs_p_group(ax, net, m):
    """Plot zone length vs probability of group

    Parameters
    ----------
    ax : axis object
    net : dict
        Dictionary with network configuration parameters
    m : int
        Group size
    """
    p = net["connection_prob"]
    p_e = net["ensemble_participation_prob"]
    N = net["ensemble_size"]
    L = net["dendrite_length"]
    Z = np.linspace(1, 100, 11) * 1e-6
    t_pop_pre = net["population_size"]
    sigma = L / int(p * t_pop_pre)
    R = net["background_rate"]
    D = net["input_duration"]

    # True connected group
    # c = 0.7384796
    c = 2.56141683
    p_total_exp_1 = []
    for z in Z:
        egrp_exp_Ei = p * N * z / L
        p_grp_exp = np.power(1 - np.exp(-egrp_exp_Ei), m)
        p_total_exp_1.append(1 - np.power((1 - p_grp_exp), c * int(L / z)))
    ax.plot(Z * 1e6, p_total_exp_1, "-o", c="b", label="true connected", alpha=0.7)

    # True active group
    # c = 0.74330444
    c = 2.57658637
    p_total_exp_1 = []
    for z in Z:
        egrp_exp_Ei = p * p_e * N * z / L
        p_grp_exp = np.power(1 - np.exp(-egrp_exp_Ei), m)
        p_total_exp_1.append(1 - np.power((1 - p_grp_exp), c * int(L / z)))
    ax.plot(Z * 1e6, p_total_exp_1, ":o", c="b", label="true active", alpha=0.7)

    # Degenerate connected group
    # c = 0.66733867
    c = 2.62999896
    p_total_exp_1 = []
    for z in Z:
        ealt_all_exp_EM = p * N * m * z / L
        p_alt_exp = 1 - np.sum(
            [poisson.pmf(k=j, mu=ealt_all_exp_EM) for j in range(0, m)]
        )
        p_total_exp_1.append(1 - np.power((1 - p_alt_exp), c * int(L / z)))
    ax.plot(
        Z * 1e6, p_total_exp_1, "-s", c="c", label="degenerate connected", alpha=0.7
    )

    # Degenerate active group
    # c = 0.67088243
    c = 2.63552124
    p_total_exp_1 = []
    for z in Z:
        ealt_all_exp_EM = p * p_e * N * m * z / L
        p_alt_exp = 1 - np.sum(
            [poisson.pmf(k=j, mu=ealt_all_exp_EM) for j in range(0, m)]
        )
        p_total_exp_1.append(1 - np.power((1 - p_alt_exp), c * int(L / z)))
    ax.plot(Z * 1e6, p_total_exp_1, ":s", c="c", label="degenerate active", alpha=0.7)

    # Noise group
    # c = 0.39755312
    c = 2.30756123
    p_total_exp_1 = []
    for z in Z:
        sigma = L / int(p * t_pop_pre)
        ENoise_exp = R * D * z * ((1 / sigma) - (p * N * m / L))
        p_noise_exp = 1 - np.sum([poisson.pmf(k=j, mu=ENoise_exp) for j in range(0, m)])
        p_total_exp_1.append(1 - np.power((1 - p_noise_exp), c * int(L / z)))
    ax.plot(Z * 1e6, p_total_exp_1, ":+", c="r", label="noise", alpha=0.7)

    # Any group
    # c = 0.32892748
    c = 2.10499296
    p_total_exp_1 = []
    for z in Z:
        sigma = L / int(p * t_pop_pre)
        EAll_exp = R * D * z * ((1 / sigma) - (p * N * m / L)) + (
            p * p_e * N * m * z / L
        )
        p_any_exp = 1 - np.sum([poisson.pmf(k=j, mu=EAll_exp) for j in range(0, m)])
        p_total_exp_1.append(1 - np.power((1 - p_any_exp), c * int(L / z)))
    ax.plot(Z * 1e6, p_total_exp_1, ":x", c="k", label="any", alpha=0.7)

    # # True connected group
    # p_total_exp_fit = []
    # for z in Z:
    #     q = N / t_pop_pre
    #     x = q * L / sigma
    #     y = L / z
    #     A = np.zeros((m + 1, m + 1))
    #     for j in range(m):
    #         A[j, j] = -(m - j) * x - j * y
    #         A[j + 1, j] = (m - j) * x
    #     for j in range(m - 1):
    #         A[j, j + 1] = (j + 1) * y
    #     B = expm(A)
    #     p_total_exp_fit.append(B[m, 0])
    # ax.plot(Z * 1e6, p_total_exp_fit, "-o", c="b", label="true connected")
    #
    # # True active group
    # p_total_exp_fit = []
    # for z in Z:
    #     q = p_e * N / t_pop_pre
    #     x = q * L / sigma
    #     y = L / z
    #     A = np.zeros((m + 1, m + 1))
    #     for j in range(m):
    #         A[j, j] = -(m - j) * x - j * y
    #         A[j + 1, j] = (m - j) * x
    #     for j in range(m - 1):
    #         A[j, j + 1] = (j + 1) * y
    #     B = expm(A)
    #     p_total_exp_fit.append(B[m, 0])
    # ax.plot(Z * 1e6, p_total_exp_fit, ":o", c="b", label="true active")
    #
    # # Degenerate connected group
    # p_total_exp_fit = []
    # for z in Z:
    #     q = (m * N) / t_pop_pre
    #     x = q * L / sigma
    #     y = L / z
    #     A = np.zeros((m + 1, m + 1))
    #     for j in range(m):
    #         A[j, j] = -x - j * y
    #         A[j + 1, j] = x
    #     for j in range(m - 1):
    #         A[j, j + 1] = (j + 1) * y
    #     B = expm(A)
    #     p_total_exp_fit.append(B[m, 0])
    # ax.plot(Z * 1e6, p_total_exp_fit, "-s", c="c", label="degenerate connected")
    #
    # # Degenerate active group
    # p_total_exp_fit = []
    # for z in Z:
    #     q = (p_e * m * N) / t_pop_pre
    #     x = q * L / sigma
    #     y = L / z
    #     A = np.zeros((m + 1, m + 1))
    #     for j in range(m):
    #         A[j, j] = -x - j * y
    #         A[j + 1, j] = x
    #     for j in range(m - 1):
    #         A[j, j + 1] = (j + 1) * y
    #     B = expm(A)
    #     p_total_exp_fit.append(B[m, 0])
    # ax.plot(Z * 1e6, p_total_exp_fit, ":s", c="c", label="degenerate active")
    #
    # # Noise group
    # p_total_exp_fit = []
    # for z in Z:
    #     q = ((t_pop_pre - m * N) * (1 - np.exp(-R * D))) / t_pop_pre
    #     x = q * L / sigma
    #     y = L / z
    #     A = np.zeros((m + 1, m + 1))
    #     for j in range(m):
    #         A[j, j] = -x - j * y
    #         A[j + 1, j] = x
    #     for j in range(m - 1):
    #         A[j, j + 1] = (j + 1) * y
    #     B = expm(A)
    #     p_total_exp_fit.append(B[m, 0])
    # ax.plot(Z * 1e6, p_total_exp_fit, ":+", c="r", label="noise")
    #
    # # Any group
    # p_total_exp_fit = []
    # for z in Z:
    #     q = (p_e * m * N + (t_pop_pre - m * N) * (1 - np.exp(-R * D))) / t_pop_pre
    #     x = q * L / sigma
    #     y = L / z
    #     A = np.zeros((m + 1, m + 1))
    #     for j in range(m):
    #         A[j, j] = -x - j * y
    #         A[j + 1, j] = x
    #     for j in range(m - 1):
    #         A[j, j + 1] = (j + 1) * y
    #     B = expm(A)
    #     p_total_exp_fit.append(B[m, 0])
    # ax.plot(Z * 1e6, p_total_exp_fit, ":x", c="k", label="any")

    sns.despine()
    ax.set_xlabel("Zone length (\u03bcm)")
    ax.set_yscale("log", nonpositive="mask")
    ax.set_ylim([1e-5, 1e0])
    ax.set_ylabel("p_group")
    # ax.legend(frameon=False)


def plot_ensemble_size_vs_p_group(ax, net, m):
    """Plot ensemble size vs probability of group

    Parameters
    ----------
    ax : axis object
    net : dict
        Dictionary with network configuration parameters
    m : int
        Group size
    """
    p = net["connection_prob"]
    p_e = net["ensemble_participation_prob"]
    N = np.logspace(1, 4, num=7)
    L = net["dendrite_length"]
    Z = net["zone_length"]
    t_pop_pre = net["population_size"]
    sigma = L / int(p * t_pop_pre)
    R = net["background_rate"]
    D = net["input_duration"]

    # True connected group
    # c = 0.7384796
    c = 2.56141683
    p_total_exp_1 = []
    for n in N:
        egrp_exp_Ei = p * n * Z / L
        p_grp_exp = np.power(1 - np.exp(-egrp_exp_Ei), m)
        p_total_exp_1.append(1 - np.power((1 - p_grp_exp), c * int(L / Z)))
    ax.plot(N, p_total_exp_1, "-o", c="b", label="true connected", alpha=0.7)

    # True active group
    # c = 0.74330444
    c = 2.57658637
    p_total_exp_1 = []
    for n in N:
        egrp_exp_Ei = p * p_e * n * Z / L
        p_grp_exp = np.power(1 - np.exp(-egrp_exp_Ei), m)
        p_total_exp_1.append(1 - np.power((1 - p_grp_exp), c * int(L / Z)))
    ax.plot(N, p_total_exp_1, ":o", c="b", label="true active", alpha=0.7)

    # Degenerate connected group
    # c = 0.66733867
    c = 2.62999896
    p_total_exp_1 = []
    for n in N:
        ealt_all_exp_EM = p * n * m * Z / L
        p_alt_exp = 1 - np.sum(
            [poisson.pmf(k=j, mu=ealt_all_exp_EM) for j in range(0, m)]
        )
        p_total_exp_1.append(1 - np.power((1 - p_alt_exp), c * int(L / Z)))
    ax.plot(N, p_total_exp_1, "-s", c="c", label="degenerate connected", alpha=0.7)

    # Degenerate active group
    # c = 0.67088243
    c = 2.63552124
    p_total_exp_1 = []
    for n in N:
        ealt_all_exp_EM = p * p_e * n * m * Z / L
        p_alt_exp = 1 - np.sum(
            [poisson.pmf(k=j, mu=ealt_all_exp_EM) for j in range(0, m)]
        )
        p_total_exp_1.append(1 - np.power((1 - p_alt_exp), c * int(L / Z)))
    ax.plot(N, p_total_exp_1, ":s", c="c", label="degenerate active", alpha=0.7)

    # Noise group
    # c = 0.39755312
    c = 2.30756123
    p_total_exp_1 = []
    for n in N:
        sigma = L / int(p * t_pop_pre)
        ENoise_exp = R * D * Z * ((1 / sigma) - (p * n * m / L))
        p_noise_exp = 1 - np.sum([poisson.pmf(k=j, mu=ENoise_exp) for j in range(0, m)])
        p_total_exp_1.append(1 - np.power((1 - p_noise_exp), c * int(L / Z)))
    ax.plot(N, p_total_exp_1, ":+", c="r", label="noise", alpha=0.7)

    # Any group
    # c = 0.32892748
    c = 2.10499296
    p_total_exp_1 = []
    for n in N:
        sigma = L / int(p * t_pop_pre)
        EAll_exp = R * D * Z * ((1 / sigma) - (p * n * m / L)) + (
            p * p_e * n * m * Z / L
        )
        p_any_exp = 1 - np.sum([poisson.pmf(k=j, mu=EAll_exp) for j in range(0, m)])
        p_total_exp_1.append(1 - np.power((1 - p_any_exp), c * int(L / Z)))
    ax.plot(N, p_total_exp_1, ":x", c="k", label="any", alpha=0.7)

    # # True connected group
    # p_total_exp_fit = []
    # for n in N:
    #     q = n / t_pop_pre
    #     x = q * L / sigma
    #     y = L / Z
    #     A = np.zeros((m + 1, m + 1))
    #     for j in range(m):
    #         A[j, j] = -(m - j) * x - j * y
    #         A[j + 1, j] = (m - j) * x
    #     for j in range(m - 1):
    #         A[j, j + 1] = (j + 1) * y
    #     B = expm(A)
    #     p_total_exp_fit.append(B[m, 0])
    # ax.plot(N, p_total_exp_fit, "-o", c="b", label="true connected")
    #
    # # True active group
    # p_total_exp_fit = []
    # for n in N:
    #     q = p_e * n / t_pop_pre
    #     x = q * L / sigma
    #     y = L / Z
    #     A = np.zeros((m + 1, m + 1))
    #     for j in range(m):
    #         A[j, j] = -(m - j) * x - j * y
    #         A[j + 1, j] = (m - j) * x
    #     for j in range(m - 1):
    #         A[j, j + 1] = (j + 1) * y
    #     B = expm(A)
    #     p_total_exp_fit.append(B[m, 0])
    # ax.plot(N, p_total_exp_fit, ":o", c="b", label="true active")
    #
    # # Degenerate connected group
    # p_total_exp_fit = []
    # for n in N:
    #     q = (m * n) / t_pop_pre
    #     x = q * L / sigma
    #     y = L / Z
    #     A = np.zeros((m + 1, m + 1))
    #     for j in range(m):
    #         A[j, j] = -x - j * y
    #         A[j + 1, j] = x
    #     for j in range(m - 1):
    #         A[j, j + 1] = (j + 1) * y
    #     B = expm(A)
    #     p_total_exp_fit.append(B[m, 0])
    # ax.plot(N, p_total_exp_fit, "-s", c="c", label="degenerate connected")
    #
    # # Degenerate active group
    # p_total_exp_fit = []
    # for n in N:
    #     q = (p_e * m * n) / t_pop_pre
    #     x = q * L / sigma
    #     y = L / Z
    #     A = np.zeros((m + 1, m + 1))
    #     for j in range(m):
    #         A[j, j] = -x - j * y
    #         A[j + 1, j] = x
    #     for j in range(m - 1):
    #         A[j, j + 1] = (j + 1) * y
    #     B = expm(A)
    #     p_total_exp_fit.append(B[m, 0])
    # ax.plot(N, p_total_exp_fit, ":s", c="c", label="degenerate active")
    #
    # # Noise group
    # p_total_exp_fit = []
    # for n in N:
    #     q = ((t_pop_pre - m * n) * (1 - np.exp(-R * D))) / t_pop_pre
    #     x = q * L / sigma
    #     y = L / Z
    #     A = np.zeros((m + 1, m + 1))
    #     for j in range(m):
    #         A[j, j] = -x - j * y
    #         A[j + 1, j] = x
    #     for j in range(m - 1):
    #         A[j, j + 1] = (j + 1) * y
    #     B = expm(A)
    #     p_total_exp_fit.append(B[m, 0])
    # ax.plot(N, p_total_exp_fit, ":+", c="r", label="noise")
    #
    # # Any group
    # p_total_exp_fit = []
    # for n in N:
    #     q = (p_e * m * n + (t_pop_pre - m * n) * (1 - np.exp(-R * D))) / t_pop_pre
    #     x = q * L / sigma
    #     y = L / Z
    #     A = np.zeros((m + 1, m + 1))
    #     for j in range(m):
    #         A[j, j] = -x - j * y
    #         A[j + 1, j] = x
    #     for j in range(m - 1):
    #         A[j, j + 1] = (j + 1) * y
    #     B = expm(A)
    #     p_total_exp_fit.append(B[m, 0])
    # ax.plot(N, p_total_exp_fit, ":x", c="k", label="any")

    sns.despine()
    ax.set_xlabel("Ensemble size N")
    ax.set_xscale("log", nonpositive="mask")
    ax.set_yscale("log", nonpositive="mask")
    ax.set_ylim([1e-5, 1e0])
    ax.set_ylabel("p_group")
    # ax.legend(frameon=False)


def plot_p_noise_group_vs_noise(ax, net, m):
    """Plot probability of noise group as a function of the rate of background activity and input duration

    Parameters
    ----------
    ax : axis object
    net : dict
        Dictionary with network configuration parameters
    m : int
        Group size
    """
    p = net["connection_prob"]
    N = net["ensemble_size"]
    L = net["dendrite_length"]
    Z = net["zone_length"]
    t_pop_pre = net["population_size"]
    sigma = L / int(p * t_pop_pre)
    R_values = np.logspace(-3, 1, num=9)
    D_values = np.logspace(-3, 1, num=9)
    c = 2.30756123

    # Noise group
    p_total_exp_fit = np.zeros((len(D_values), len(R_values)))
    for d, D in enumerate(D_values):
        for r, R in enumerate(R_values):
            ENoise_exp = R * D * Z * ((1 / sigma) - (p * N * m / L))
            p_noise_exp = 1 - np.sum(
                [poisson.pmf(k=j, mu=ENoise_exp) for j in range(0, m)]
            )
            p_total_exp_fit[r, d] = 1 - np.power((1 - p_noise_exp), c * int(L / Z))
            # q = ((t_pop_pre - m * N) * (1 - np.exp(-R * D))) / t_pop_pre
            # x = q * L / sigma
            # y = L / Z
            # A = np.zeros((m + 1, m + 1))
            # for j in range(m):
            #     A[j, j] = -x - j * y
            #     A[j + 1, j] = x
            # for j in range(m - 1):
            #     A[j, j + 1] = (j + 1) * y
            # B = expm(A)
            # p_total_exp_fit[r, d] = B[m, 0]

    # Clip values at 1e-6
    p_total_exp_fit = np.clip(p_total_exp_fit, a_min=1e-6, a_max=None)
    ax = sns.heatmap(
        p_total_exp_fit,
        ax=ax,
        norm=LogNorm(1e-6, 1e0),
        xticklabels=R_values,
        yticklabels=D_values,
        cmap=cm.viridis,
        cbar_kws={
            "extend": "min",
            "label": "p_noise_group",
            "ticks": [1e-6, 1e-3, 1e0],
            "format": mticker.FixedFormatter(["<=$10^{-6}$", "$10^{-3}$", "$10^{0}$"]),
        },
        square=True,
    )

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


def plot_p_any_group_vs_noise(ax, net, m):
    """Plot probability of any group as a function of the rate of background activity and input duration

    Parameters
    ----------
    ax : axis object
    net : dict
        Dictionary with network configuration parameters
    m : int
        Group size
    """
    p = net["connection_prob"]
    p_e = net["ensemble_participation_prob"]
    N = net["ensemble_size"]
    L = net["dendrite_length"]
    Z = net["zone_length"]
    t_pop_pre = net["population_size"]
    sigma = L / int(p * t_pop_pre)
    R_values = np.logspace(-3, 1, num=9)
    D_values = np.logspace(-3, 1, num=9)

    # Any group
    p_total_exp_fit = np.zeros((len(D_values), len(R_values)))
    c = 2.10499296
    for d, D in enumerate(D_values):
        for r, R in enumerate(R_values):
            EAll_exp = R * D * Z * ((1 / sigma) - (p * N * m / L)) + (
                p * p_e * N * m * Z / L
            )
            p_any_exp = 1 - np.sum([poisson.pmf(k=j, mu=EAll_exp) for j in range(0, m)])
            p_total_exp_fit[r, d] = 1 - np.power((1 - p_any_exp), c * int(L / Z))
            # q = (p_e * m * N + (t_pop_pre - m * N) * (1 - np.exp(-R * D))) / t_pop_pre
            # x = q * L / sigma
            # y = L / Z
            # A = np.zeros((m + 1, m + 1))
            # for j in range(m):
            #     A[j, j] = -x - j * y
            #     A[j + 1, j] = x
            # for j in range(m - 1):
            #     A[j, j + 1] = (j + 1) * y
            # B = expm(A)
            # p_total_exp_fit[r, d] = B[m, 0]
    ax = sns.heatmap(
        p_total_exp_fit,
        ax=ax,
        norm=LogNorm(1e-6, 1e0),
        xticklabels=R_values,
        yticklabels=D_values,
        cmap=cm.viridis,
        cbar_kws={
            "label": "p_any_group",
            "ticks": [1e-6, 1e-3, 1e0],
        },
        square=True,
    )

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


def plot_legend(ax, bx):
    """Plot legend

    Parameters
    ----------
    ax : axis object
    """
    ax.text(
        0.4,
        0.9,
        "Legend for A, B, C, D",
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )
    handles, labels = bx.get_legend_handles_labels()
    ax.legend(handles=handles, labels=labels, loc="center", frameon=False)
    ax.axis("off")


def main():
    figure = plt.figure(figsize=(10, 10), constrained_layout=True)
    spec = figure.add_gridspec(3, 8)
    axa = figure.add_subplot(spec[0, 0:4])
    axb = figure.add_subplot(spec[0, 4:8])
    # axc = figure.add_subplot(spec[0, 6:8])
    axd = figure.add_subplot(spec[1, 0:4])
    axe = figure.add_subplot(spec[1, 4:8])
    # axf = figure.add_subplot(spec[1, 6:8])
    axg = figure.add_subplot(spec[2, 0:4])
    axh = figure.add_subplot(spec[2, 4:8])

    plot_conn_prob_vs_p_group(axa, networks[4], 5)
    plot_ensemble_participation_prob_vs_p_group(axb, networks[4], 5)
    plot_zone_length_vs_p_group(axd, networks[4], 5)
    plot_ensemble_size_vs_p_group(axe, networks[4], 5)
    plot_p_noise_group_vs_noise(axg, networks[4], 5)
    plot_p_any_group_vs_noise(axh, networks[4], 5)
    # plot_legend(axc, axa)

    plt.figtext(0.01, 0.98, "A", fontsize=12, weight="bold")
    plt.figtext(0.52, 0.98, "B", fontsize=12, weight="bold")
    plt.figtext(0.01, 0.65, "C", fontsize=12, weight="bold")
    plt.figtext(0.52, 0.65, "D", fontsize=12, weight="bold")
    plt.figtext(0.06, 0.33, "E", fontsize=12, weight="bold")
    plt.figtext(0.56, 0.33, "F", fontsize=12, weight="bold")

    plt.savefig("Figure-2-param-screen.png", bbox_inches="tight")
    plt.savefig("Figure-2-param-screen.svg", bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
