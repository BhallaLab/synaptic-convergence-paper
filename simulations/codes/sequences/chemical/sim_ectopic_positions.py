#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : sim_ectopic_positions.py
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
    "numDesiredSeqPermutations": 24,  # number of sequence permutations
    "numExtra": 8,  # number of extra inputs
}


def buildModel(chemName, stim, extraPlace, extraTime):
    """Build model

    Parameters
    ----------
    chemName : str
        Name of model
    stim : array
        Time order of arrival of inputs
    extraPlace : int
        Voxel where additional input arrives relative to the start position of the first stim input
    extraTime : float
        Time point in units of seqDt when additional input arrives relative to the start of first stim input

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

    k = int(params["blanks"] + extraPlace)
    seqExtra[k].nInit = params["preStim"] + extraTime * params["seqDt"]
    Zseq[k].nInit = 1
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


def runStimulus(stim, extraPlace, extraTime):
    """Run simulation for a stimulus pattern

    Parameters
    ----------
    stim : array
        Array of timing order of arrival of inputs
    extraPlace : int
        Voxel where additional input arrives relative to the start position of the first stim input
    extraTime : float
        Time point in multiples of seqDt when additional input arrives relative to the start of first stim input

    Returns
    -------
    auc : float
        Total area under curve
    """
    moose.Neutral("/library")
    makePassiveSoma("cell", params["dendLength"], params["dendDiameter"])
    ame.makeBis(params)
    model = buildModel("bis", stim, extraPlace, extraTime)
    moose.reinit()
    moose.start(
        params["preStim"] + params["postStim"] + params["seqDt"] * params["numInputs"]
    )
    auc = extractAuc("/model/graphs/plot0")
    moose.delete(model)
    moose.delete("/library")
    return auc


def calcSelectivity(stimArray, extraPlaceList, extraTime):
    """Run different stimuli and calculate selectivity

    Parameters
    ----------
    stimArray : ndarray
        List of stimuli where each stimulus contains the time order of arrival of inputs
    extraPlaceList : array
        List of positions where additional inputs arrive
    extraTime : float
        Time point in multiples of seqDt when additional input arrives relative to the start of first stim input

    Returns
    -------
    selectivity : array
        Array of selectivity values where selectivity is measured as (response(sequence) - mean_response(all_sequences)) / max_response
    """
    selectivity = []
    for extraPlace in extraPlaceList:
        aucs = []
        for stim in stimArray:
            aucs.append(runStimulus(stim, extraPlace, extraTime))
        selectivity.append((aucs[0] - np.mean(aucs)) / max(aucs))
    return selectivity


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

    stride = int(params["dist"] / params["numInputs"])
    extraPlaceList = np.arange(
        -params["numExtra"], (params["numInputs"] - 1) * stride + 1 + params["numExtra"]
    )
    extraPlaceList = np.delete(
        extraPlaceList,
        params["numExtra"]
        + np.arange(0, (params["numInputs"] - 1) * stride + 1, stride),
    )
    extraTimeList = [0, int(params["numInputs"] / 2), params["numInputs"] - 1]
    timeLabels = ["t_start", "t_mid", "t_end"]
    extraTimeDict = dict(zip(timeLabels, extraTimeList))
    outputs = {}
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))
    outputs["positions"] = extraPlaceList * params["diffusionLength"] * 1e6
    outputs["zone"] = (0, (params["numInputs"] - 1) * stride)
    reference_selectivity = calcSelectivity(stimArray, [0], 10000)
    outputs["reference"] = reference_selectivity
    ax.axhline(reference_selectivity, color="k", label="reference")
    for label, et in extraTimeDict.items():
        selectivity = calcSelectivity(stimArray, extraPlaceList, et)
        outputs[label] = selectivity
        ax.plot(
            extraPlaceList * params["diffusionLength"] * 1e6,
            selectivity,
            "-",
            label=label,
        )

    np.save(DATA_PATH + "ectopic_positions_selectivity.npy", outputs)
    ax.axvspan(
        0,
        (params["numInputs"] - 1) * stride * params["diffusionLength"] * 1e6,
        color="c",
        alpha=0.5,
        lw=0,
    )
    ax.set_xlabel("Position of extra input (\u03bcm)")
    ax.set_ylabel("Selectivity")
    plt.legend()
    plt.savefig("ectopic_positions_selectivity.png")
    plt.show()


if __name__ == "__main__":
    main()
