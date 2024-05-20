#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : sim_selectivty_seq_length.py
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
import math

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
    "numDesiredStimPermutations": 24,  # number of stim permutations
    "numExtra": 8,  # number of extra inputs
    "seqDx": 3e-6,  # Spacing between inputs
}


def buildModel(chemName, stim, extraPlace, extraTime, doExtra2=False):
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
    doExtra2 : bool
        Give second extra input

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
    seqphase.nInit = 10000  # Some large value greater than the simulation duration
    seqExtra.nInit = 10000  # Some large value greater than the simulation duration
    Zseq.nInit = 0
    for j in range(params["numInputs"]):
        k = int(params["blanks"] + j * params["stride"])
        Zseq[k].nInit = 1
        seqphase[k].nInit = params["preStim"] + stim[j] * params["seqDt"]

    k = int(params["blanks"] + extraPlace)
    seqExtra[k].nInit = params["preStim"] + extraTime * params["seqDt"]
    Zseq[k].nInit = 1

    if doExtra2:
        extraPlace2 = math.ceil(((params["numInputs"] - 1) * params["stride"] + 1) / 2)
        extraTime2 = int(params["numInputs"] / 2)
        if extraPlace == extraPlace2 and extraTime == extraTime2:
            extraPlace2 -= 1
            extraTime2 -= 1
        if extraTime >= 10000:
            extraTime2 = 10000
        k = int(params["blanks"] + extraPlace2)
        seqExtra[k].nInit = params["preStim"] + extraTime2 * params["seqDt"]
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


def runStimulus(stim, extraPlace, extraTime, doExtra2):
    """Run simulation for a stimulus pattern

    Parameters
    ----------
    stim : array
        Array of timing order of arrival of inputs
    extraPlace : int
        Voxel where additional input arrives relative to the start position of the first stim input
    extraTime : float
        Time point in multiples of seqDt when additional input arrives relative to the start of first stim input
    doExtra2 : bool
        Give second extra input

    Returns
    -------
    auc : float
        Total area under curve
    """
    moose.Neutral("/library")
    makePassiveSoma("cell", params["dendLength"], params["dendDiameter"])
    ame.makeBis(params)
    model = buildModel("bis", stim, extraPlace, extraTime, doExtra2)
    moose.reinit()
    moose.start(
        params["preStim"] + params["postStim"] + params["seqDt"] * params["numInputs"]
    )
    auc = extractAuc("/model/graphs/plot0")
    moose.delete(model)
    moose.delete("/library")
    return auc


def calcSelectivity(stimArray, extraPlace, extraTime, doExtra2):
    """Run different stimuli and calculate selectivity

    Parameters
    ----------
    stimArray : ndarray
        List of stimuli where each stimulus contains the time order of arrival of inputs
    extraPlace : int
        Voxel where additional input arrives relative to the start position of the first stim input
    extraTime : float
        Time point in multiples of seqDt when additional input arrives relative to the start of first stim input
    doExtra2 : bool
        Give second extra input

    Returns
    -------
    selectivity : float
        Array of selectivity values where selectivity is measured as (response(sequence) - mean_response(all_sequences)) / max_response
    """
    aucs = []
    for stim in stimArray:
        aucs.append(runStimulus(stim, extraPlace, extraTime, doExtra2))
    selectivity = (aucs[0] - np.mean(aucs)) / max(aucs)
    return selectivity


def main():
    """Create a compartment and simulate bistable-switch system for different number of inputs"""
    global params

    if not os.path.exists(DATA_PATH):
        os.makedirs(DATA_PATH)
    outputs = {}
    outputs["reference"] = []
    outputs["1_ectopic"] = {}
    outputs["2_ectopic"] = {}
    inputNumberRange = np.arange(3, 11)
    stimAmplitude = [1.10, 1.05, 1.0, 1.0, 0.98, 0.95, 0.94, 0.94]
    for n, ampl in zip(inputNumberRange, stimAmplitude):
        params["numInputs"] = n
        params["stimAmpl"] = ampl
        outputs["1_ectopic"][n] = []
        outputs["2_ectopic"][n] = []
        params["stride"] = int(params["seqDx"] / params["diffusionLength"])
        print(params["numInputs"])
        print(params["stride"])
        perm = list(itertools.permutations(range(params["numInputs"])))
        if params["numInputs"] > 3:
            stimArray = perm[0 :: int(len(perm) / params["numDesiredStimPermutations"])]
        else:
            stimArray = perm

        print(stimArray)

        outputs["reference"].append(
            calcSelectivity(stimArray, 0, 10000, doExtra2=False)
        )
        extraPlaceList = [
            1,
            math.ceil(((params["numInputs"] - 1) * params["stride"] + 1) / 2),
            (params["numInputs"] - 1) * params["stride"] - 1,
        ]
        print(extraPlaceList)
        extraTimeList = [0, int(params["numInputs"] / 2), params["numInputs"]]
        print(extraTimeList)
        for ep in extraPlaceList:
            for et in extraTimeList:
                outputs["1_ectopic"][n].append(
                    calcSelectivity(stimArray, ep, et, doExtra2=False)
                )
                outputs["2_ectopic"][n].append(
                    calcSelectivity(stimArray, ep, et, doExtra2=True)
                )

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))
    ax.plot(inputNumberRange, outputs["reference"], color="k", label="reference")
    ax.plot(
        inputNumberRange,
        [np.mean(outputs["1_ectopic"][i]) for i in outputs["1_ectopic"].keys()],
        "-",
        label="1_ectopic",
    )
    ax.plot(
        inputNumberRange,
        [np.mean(outputs["2_ectopic"][i]) for i in outputs["2_ectopic"].keys()],
        "-",
        label="2_ectopic",
    )
    ax.legend(frameon=False)
    ax.set_xlabel("Sequence Length")
    ax.set_ylabel("Selectivity")

    np.save(DATA_PATH + "seq_lengths_selectivity.npy", outputs)
    ax.set_xlabel("Sequence length (\u03bcm)")
    ax.set_ylabel("Selectivity")
    plt.legend()
    plt.savefig("seq_lengths_selectivity.png")
    plt.show()


if __name__ == "__main__":
    main()
