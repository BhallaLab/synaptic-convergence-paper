#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : sim_seq_vs_scrambled.py
# Author            : Bhanu Priya S <bhanups@ncbs.res.in>
# Date              : 07.02.2023
# Last Modified Date: 09.11.2023
# Last Modified By  : Bhanu Priya S <bhanups@ncbs.res.in>

import numpy as np
import matplotlib.pyplot as plt
import moose
import abstrModelEqns10 as ame
import rdesigneur as rd
import os
import itertools
import csv
from scipy.stats import linregress

DATA_PATH = "./data/"  # Path to store data

# All params in SI units
params = {
    "dendDiameter": 10e-6,  # Diameter of section of dendrite in model, voxel length
    "dendLength": 100e-6,  # Length of section of dendrite in model
    "diffusionLength": 1e-6,  # Spatial discretiztion for diffusion
    "chemPlotDt": 0.1,  # Timestep for storing and plotting chemical values
    "chemDt": 0.02,  # Timestep for chemical computations
    "diffDt": 0.002,  # Timestep for diffusion
    "preStim": 10.0,  # pre-stim duration
    "postStim": 80.0,  # post-stim duration
    "blanks": 40,  # number of blank voxels
    "numInputs": 5,  # number of inputs in stimulus
    "dist": 15.0,  # Zone length over which inputs arrive
    "seqDt": 2.0,  # Time interval between subsequent inputs
    "stimAmpl": 1,  # Scaling constant for stimulus
    "numDesiredSeqPermutations": 120,  # number of sequence permutations
}


def buildModel(chemName, stim):
    """Build model

    Parameters
    ----------
    chemName : str
        Name of model
    stim : array
        Time order of arrival of inputs

    Returns
    -------
    rdes.model : MOOSE model
        MOOSE model
    """
    rdes = rd.rdesigneur(
        useGssa=False,
        turnOffElec=True,
        verbose=False,
        chemPlotDt=params["chemPlotDt"],
        chemDt=params["chemDt"],
        diffDt=params["diffDt"],
        diffusionLength=params["diffusionLength"],
        cellProto=[["cell", "soma"]],
        chemProto=[["dend", chemName]],
        chemDistrib=[["dend", "soma", "install", "1"]],
        plotList=[["soma", "1", "dend" + "/A", "n", "# of A"]],
    )
    moose.element("/library/bis").name = "chem"
    moose.element("/library/chem/bis").name = "dend"
    rdes.buildModel()
    Aseq = moose.vec("/model/chem/dend/A")
    CaSeq = moose.vec("/model/chem/dend/Ca")
    Zseq = moose.vec("/model/chem/dend/Z")
    seqphase = moose.vec("/model/chem/dend/phase")
    seqExtra = moose.vec("/model/chem/dend/extra")
    moose.vec("/model/chem/dend/ampl").nInit = params["stimAmpl"]
    stride = int(params["dist"] / params["numInputs"])
    seqphase.nInit = 10000  # Some large value greater than the simulation duration
    seqExtra.nInit = 10000  # Some large value greater than the simulation duration
    Zseq.nInit = 0
    for j in range(params["numInputs"]):
        k = int(params["blanks"] + j * stride)
        Zseq[k].nInit = 1
        seqphase[k].nInit = params["preStim"] + stim[j] * params["seqDt"]
    return rdes.model


def extractAuc(tabname):
    """Extract area under curve

    Parameters
    ----------
    tabname : str
        Table name

    Returns
    -------
    auc : float
        Total area under curve
    """
    tabs = moose.vec(tabname)
    auc = 0.0
    for i in range(len(tabs)):
        data = tabs[i].vector
        auc += np.sum(data)
    return auc


def makePassiveSoma(name, length, diameter):
    """Create a passive soma with a dendrite

    Parameters
    ----------
    name : str
        Compartment name
    length : float
        Length of the compartment
    diameter : float
        Diameter of the compartment

    Returns
    -------
    elecid : Neuron
        Neuron object with the dendrite
    """
    elecid = moose.Neuron("/library/" + name)
    dend = moose.Compartment(elecid.path + "/soma")
    dend.diameter = diameter
    dend.length = length
    dend.x = length
    return elecid


def runStimulus(stim):
    """Run simulation for a stimulus pattern

    Parameters
    ----------
    stim : array
        Array of timing order of arrival of inputs

    Returns
    -------
    auc : float
        Total area under curve
    """

    moose.Neutral("/library")
    makePassiveSoma("cell", params["dendLength"], params["dendDiameter"])
    ame.makeBis(params)
    model = buildModel("bis", stim)
    moose.reinit()
    moose.start(
        params["preStim"] + params["postStim"] + params["seqDt"] * params["numInputs"]
    )
    auc = extractAuc("/model/graphs/plot0")
    moose.delete(model)
    moose.delete("/library")
    return auc


def main():
    """Create a compartment and simulate bistable-switch system"""
    global params

    if not os.path.exists(DATA_PATH):
        os.makedirs(DATA_PATH)
    perm = list(itertools.permutations(range(params["numInputs"])))
    if params["numInputs"] > 3:
        stimArray = perm[0 :: int(len(perm) / params["numDesiredSeqPermutations"])]
    else:
        stimArray = perm

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))
    aucs = []
    Q_scores = []
    for stim in stimArray:
        aucs.append(runStimulus(stim))
        slope, intercept, r, p, std_err = linregress(np.arange(len(stim)), stim)
        Q_scores.append(slope * r * r)
    ax.scatter(Q_scores[1:], aucs[1:], label="other patterns")
    ax.scatter(
        Q_scores[0],
        aucs[0],
        marker="x",
        color="r",
        label="ordered sequence",
    )
    ax.set_xlabel("Directionality (Q score)")
    ax.set_ylabel("A total (a.u.)")
    ax.set_xlim([-1.05, 1.05])
    ax.set_ylim(ymin=0)
    ax.legend(frameon=False)
    outputs = {}
    outputs["stim_pattern"] = stimArray
    outputs["Q_score"] = Q_scores
    outputs["A_total"] = aucs
    np.save(DATA_PATH + "responses_to_stim_patterns.npy", outputs)
    plt.savefig("responses_to_stim_patterns.png")
    plt.show()


if __name__ == "__main__":
    main()
