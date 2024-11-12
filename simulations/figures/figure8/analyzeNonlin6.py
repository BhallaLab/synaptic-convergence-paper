import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt

AMPAR_gmax = 0.157  # nS
NMDAR_gmax = 0.063  # nS


def extractStats(fname):
    # Read from HDF5
    eSettleIdx = 20000
    cSettleIdx = 400
    with pd.HDFStore(fname, mode="r") as store:
        # Read the DataFrame
        spinedf = store["spineCa"].values * 1000
        camdf = store["dendCaM"].values * 1000
        mapkdf = store["dendMAPK"].values * 1e6
        Vm = store["somaVm"].values * 1000
        branchVm = store["branchVm"].values * 1000
        metadata = store.get_storer("dendCaM").attrs.metadata
        chemDt = metadata["chemDt"]
        runtime = metadata["runtime"]

        spineBaseline = np.mean(spinedf[cSettleIdx // 2 : cSettleIdx])
        camBaseline = np.mean(camdf[cSettleIdx // 2 : cSettleIdx])
        mapkBaseline = np.mean(mapkdf[cSettleIdx // 2 : cSettleIdx])

        elecDt = runtime / (len(Vm) - 1)
        # print( len( Vm ), runtime, elecDt, metadata['elecPlotDt'], chemDt, spinedf.shape, camdf.shape, mapkdf.shape )
        # print( spineBaseline, camBaseline, mapkBaseline )

        return [
            np.max(Vm[eSettleIdx:] + 72),
            np.max(branchVm[eSettleIdx:] + 72),
            np.mean(spinedf[cSettleIdx:,]) - spineBaseline,
            np.max(camdf[cSettleIdx:,]) - camBaseline,
            np.mean(mapkdf[cSettleIdx:,]) - mapkBaseline,
        ]


def doPlot(ax, xtab, ytab, panel, xlabel, ylabel):
    print("Plotting: ", panel, ylabel, xlabel)
    ax.plot(xtab, ytab)
    ax.tick_params(axis="both", which="major", labelsize=14)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    # ax.set_title( title )
    ax.set_ylim(0, max(ytab) * 1.1)
    # ax.legend( fontsize = 14, frameon = False )
    ax.text(-0.15, 1.05, panel, fontsize=14, weight="bold", transform=ax.transAxes)
    # ax.set_ylim( -75, -55 )


def main():
    vmIdx = 0
    camIdx = 3
    mapkIdx = 4
    fig = plt.figure(figsize=(10, 7))
    gs = fig.add_gridspec(2, 3)

    ################################################
    xtab = np.array([5, 10, 15, 20, 25, 30, 35, 40])
    ytab = []
    for wt in xtab:
        ret = extractStats("out7_frq1.0_Dt0.01_Dx6_N7_seq_p0_w{}.0.h5".format(wt))
        ytab.append(ret[vmIdx])

    ax = fig.add_subplot(gs[1, 0])
    doPlot(
        ax,
        xtab * AMPAR_gmax,
        ytab,
        "C",
        "$g_{max\_AMPAR}$ (nS)",
        "Vm (mV)",
    )
    ax2 = ax.twiny()
    doPlot(
        ax2,
        xtab * NMDAR_gmax,
        ytab,
        "C",
        "$g_{max\_NMDAR}$ (nS)",
        "Vm (mV)",
    )
    ax2.spines["top"].set_visible(True)
    ################################################

    xtab = [1, 2, 3, 4, 5, 6, 8, 10]
    ytab = []
    for dx in xtab:
        ret = extractStats("out7_frq1.0_Dt0.01_Dx{}_N7_seq_p0_w20.0.h5".format(dx))
        ytab.append(ret[vmIdx])
    nxtab = np.array(xtab) * 2  # Synapse Spacing is 2 microns.
    ax = fig.add_subplot(gs[0, 0])
    doPlot(ax, nxtab, ytab, "B", "Syn Spacing (microns)", "Vm (mV)")
    ################################################

    xtab = np.array([5, 10, 15, 20, 25, 30, 35, 40])
    ytab = []
    for wt in xtab:
        ret = extractStats("out7_frq20.0_Dt3.0_Dx2_N5_seq_p0_w{}.0.h5".format(wt))
        ytab.append(ret[mapkIdx])
    ax = fig.add_subplot(gs[1, 2])
    doPlot(
        ax,
        xtab * AMPAR_gmax,
        ytab,
        "G",
        "$g_{max\_AMPAR}$ (nS)",
        "MAPK-p (AUC)",
    )
    ax2 = ax.twiny()
    doPlot(
        ax2,
        xtab * NMDAR_gmax,
        ytab,
        "G",
        "$g_{max\_NMDAR}$ (nS)",
        "MAPK-p (AUC)",
    )
    ax2.spines["top"].set_visible(True)
    ################################################

    xtab = [1, 2, 3, 4, 5, 6, 8, 10]
    ytab = []
    for dx in xtab:
        ret = extractStats("out7_frq20.0_Dt3.0_Dx{}_N5_seq_p0_w20.h5".format(dx))
        ytab.append(ret[mapkIdx])
    nxtab = np.array(xtab) * 2  # Synapse Spacing is 2 microns.
    ax = fig.add_subplot(gs[0, 2])
    doPlot(ax, nxtab, ytab, "F", "Syn Spacing (microns)", "MAPK-p (AUC)")
    ################################################

    xtab = np.array([5, 10, 15, 20, 25, 30, 35, 40])
    ytab = []
    for wt in xtab:
        ret = extractStats("out7_frq10.0_Dt1.0_Dx1_N5_grp_p0_w{}.0.h5".format(wt))
        ytab.append(ret[camIdx])
    ax = fig.add_subplot(gs[1, 1])
    doPlot(
        ax,
        xtab * AMPAR_gmax,
        ytab,
        "E",
        "$g_{max\_AMPAR}$ (nS)",
        "Ca4.CaM (AUC)",
    )
    ax2 = ax.twiny()
    doPlot(
        ax2,
        xtab * NMDAR_gmax,
        ytab,
        "E",
        "$g_{max\_NMDAR}$ (nS)",
        "Ca4.CaM (AUC)",
    )
    ax2.spines["top"].set_visible(True)

    ################################################

    xtab = [1, 2, 3, 4, 5]
    ytab = []
    for dx in xtab:
        ret = extractStats("out7_frq10.0_Dt1.0_Dx{}_N5_grp_p0_w20.h5".format(dx))
        ytab.append(ret[camIdx])
    nxtab = np.array(xtab) * 2  # Synapse Spacing is 2 microns.
    ax = fig.add_subplot(gs[0, 1])
    doPlot(ax, nxtab, ytab, "D", "Syn Spacing (microns)", "Ca4.CaM (AUC)")
    ################################################

    plt.tight_layout()
    plt.savefig("nonlin6.png")
    plt.savefig("nonlin6.svg")
    plt.show()


if __name__ == "__main__":
    main()
