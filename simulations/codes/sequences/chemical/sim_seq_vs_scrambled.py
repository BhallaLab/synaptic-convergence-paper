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
import csv

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
    "numInputs": 5,  # numberL
    "dist": 15.0,  # Zone length over which inputs arrive
    "seqDt": 2.0,  # Time interval between subsequent inputs
    "stimAmpl": 1,  # Scaling constant for stimulus
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


def extractTimeSlice(tabname, time_point):
    """Extract table data for all voxels at a specific time point

    Parameters
    ----------
    tabname : str
        Table name
    time_point : float
        Time point of simulation

    Returns
    -------
    x_vector : np.ndarray
        Array of voxel numbers
    a_vector : np.ndarray
        Array of data stored in table in all voxels
    """
    tab = moose.element(tabname)
    tab_vec = moose.vec(tabname)
    time_index = int(time_point / tab.dt)
    x_vector = np.arange(len(tab_vec))
    a_vector = np.array(
        [np.array(tab_vec[i].vector)[time_index] for i in range(len(tab_vec))]
    )
    return x_vector, a_vector


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
    x_vec : np.ndarray
        Array of voxel indices
    a_vec : np.ndarray
        Array of molecule numbers
    """
    moose.Neutral("/library")
    makePassiveSoma("cell", params["dendLength"], params["dendDiameter"])
    ame.makeBis(params)
    model = buildModel(
        "bis",
        stim,
    )
    moose.reinit()
    moose.start(
        params["preStim"] + params["postStim"] + params["seqDt"] * params["numInputs"]
    )
    x_vec, a_vec = extractTimeSlice(
        "/model/graphs/plot0",
        params["preStim"] + params["seqDt"] * (params["numInputs"] - 0.8),
    )
    moose.delete(model)
    moose.delete("/library")
    return x_vec, a_vec


def main():
    """Create a compartment and simulate bistable-switch system"""
    global params

    seqs = {
        "ordered": np.array([0, 1, 2, 3, 4]),
        "scrambled": np.array([4, 1, 0, 3, 2]),
    }

    if not os.path.exists(DATA_PATH):
        os.makedirs(DATA_PATH)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))
    for label, seq in seqs.items():
        x_vec, a_vec = runStimulus(seq)
        with open(
            DATA_PATH + "%s_output.csv" % (label),
            "w",
        ) as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(a_vec)
        ax.plot(x_vec * params["diffusionLength"] * 1e6, a_vec, "-", label=label)
    ax.set_xlabel("Position (\u03bcm)")
    ax.set_ylabel("A (a.u.)")
    plt.legend()
    plt.savefig("seq_vs_scram.png")
    plt.show()


if __name__ == "__main__":
    main()
