#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File              : Ca_sim6_sample_Ca_trace.py
# Author            : Bhanu Priya S <bhanups@ncbs.res.in>
# Date              : 20.01.2023
# Last Modified Date: 23.01.2023
# Last Modified By  : Bhanu Priya S <bhanups@ncbs.res.in>

"""Ca-CaM Model simulation

This code simulates ca_m activation by Ca under two conditions:
1) Grouped inputs : wherein the stimuli arrive at nearby locations
2) Dispersed inputs : wherein the stimuli arrive at locations that 
are more spaced apart
Model is specified in Ca_CaM_11.g
"""

import os
import moose
import numpy as np
import sys
import matplotlib.pyplot as plt
import rdesigneur as rd

import csv

DATA_PATH = "./data/"  # Path to store data


params = {
    "diffusion_length": 0.5e-6,  # Diffusion characteristic length, used as voxel length too.
    "dend_diameter": 5.0e-6,  # Diameter of section of dendrite in model
    "dend_length": 100e-6,  # Length of section of dendrite in model
    "diff_const_ca": 20e-12,  # Diffusion constant of Ca
    "stim_amplitude": 0.05,  # Ca Stimulus amplitude, mM
    "pre_stim_time": 10.0,  # Time to run before turning on stimulus.
    "post_stim_time": 30.0,  # Time to run after stimulus.
    "stim_width": 0.05,  # Duration of Ca influx for each stimulus.
    "chem_model": "Ca_CaM.g",  # Chem model definition
    "seq_dt": 0.0,  # Time interval between the start of successive inputs in seq
    # Zero since the inputs are simultaneous
    "seq_dx": 2.0e-6,  # Distance between successive inputs in seq.
    "dispersion_factor": 5,  #  = Dispersed_spacing/Grouped_spacing
    "num_inputs": 5,  # Number of inputs
    "seed": 12345,  # Seed for random number generator
    "basal_ca_conc": 0.08e-3,  # Basal Ca conc in mM
}


def make_passive_soma(name, length, diameter):
    """Create a passive soma with a dendrite

    Parameters
    ----------
    name : str
        compartment name
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


def set_diff_const(element, param_name):
    """Set diffusion constant

    Parameters
    ----------
    element : str
        Molecule whose diffusion constant needs to be set
    param_name : str
        Name of parameter, i.e., diffusion constant
    """
    moose.le("/library/chem/kinetics")
    moose.le("/library/chem/kinetics/Ca")
    e = moose.element("/library/chem/kinetics/Ca/" + element)
    e.diffConst = params[param_name]


def run_stimulus(stim_q):
    """
    Initialize Ca conc stimulus and run simulation.
    Return total ca_m_Ca4 concentration

    Parameters
    ----------
    stim_q : dict
        Dictionary of {time:[position_index, Ca_input_conc]}

    Returns
    -------
    mean_conc : float
        Mean concentration of ca_m_Ca4 in the dendrite across time
    """
    Ca_input = moose.vec("/model/chem/dend/Ca/Ca_input")
    Ca = moose.vec("/model/chem/dend/Ca/Ca")
    Ca_input.concInit = 0e-6
    Ca.concInit = params["basal_ca_conc"]
    moose.reinit()

    clock = moose.element("/clock")
    for t in sorted(stim_q):
        currt = clock.currentTime
        if t > currt:
            moose.start(t - currt)
            print(f"t={t} stim={stim_q[t]}")
            for entry in stim_q[t]:
                index, conc = entry

                Ca_input[int(index)].concInit = conc
    moose.start(params["post_stim_time"])
    ca_m_Ca4 = moose.vec("/model/graphs/plot2")
    mean_conc = np.mean(np.array([i.vector for i in ca_m_Ca4]), axis=0)
    return mean_conc


def simulate_ca_cam():
    """
    Setup and simulate ca_m activation by Ca for grouped and dispersed inputs
    """
    moose.seed(int(params["seed"]))
    rdes = rd.rdesigneur(
        useGssa=False,
        turnOffElec=True,
        chemDt=0.0005,
        chemPlotDt=0.001,
        diffDt=0.0005,
        diffusionLength=params["diffusion_length"],
        cellProto=[["cell", "soma"]],
        chemProto=[[params["chem_model"], "chem"]],
        chemDistrib=[["chem", "soma", "install", "1"]],
        plotList=[
            ["soma", "1", "dend/Ca/Ca", "conc", "[dend Ca(mM)]"],
            ["soma", "1", "dend/Ca/Ca_input", "conc", "[dend Ca_input(mM)]"],
            ["soma", "1", "dend/CaM/CaM_Ca4", "conc", "[dend CaM_4(mM)]"],
        ],
    )
    # Assign parameters to the prototype model.
    set_diff_const("Ca", "diff_const_ca")
    rdes.buildModel()
    print("MODEL BUILT")

    ################################################################
    # Run and print the stimulus
    # stim_q is { time:[index, ca_conc], ... }
    ampl = params["stim_amplitude"]
    basal_ampl = 0.0  # Basal Ca Input

    # Grouped stimulus
    stim1 = {}
    for n in range(params["num_inputs"]):
        stim1[params["pre_stim_time"] + n * params["seq_dt"]] = []
        stim1[
            params["pre_stim_time"] + n * params["seq_dt"] + params["stim_width"]
        ] = []

    blank_voxels_at_end = int(
        (
            int(params["dend_length"] / params["diffusion_length"])
            - (
                int(params["seq_dx"] / params["diffusion_length"])
                * (params["num_inputs"] - 1)
                + 1
            )
        )
        / 2
    )
    for n in range(params["num_inputs"]):
        location = (
            blank_voxels_at_end + int(params["seq_dx"] / params["diffusion_length"]) * n
        )
        stim1[params["pre_stim_time"] + n * params["seq_dt"]].append((location, ampl))
        stim1[
            params["pre_stim_time"] + n * params["seq_dt"] + params["stim_width"]
        ].append((location, basal_ampl))
    grouped_tot = run_stimulus(stim1)

    with open(
        DATA_PATH
        + "sample_grouped_output_seq_dx_%.6f_seq_dt_%.3f.csv"
        % (params["seq_dx"], params["seq_dt"]),
        "w",
    ) as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(grouped_tot)

    # Dispersed stimulus
    stim2 = {}
    for n in range(params["num_inputs"]):
        stim2[params["pre_stim_time"] + n * params["seq_dt"]] = []
        stim2[
            params["pre_stim_time"] + n * params["seq_dt"] + params["stim_width"]
        ] = []

    blank_voxels_at_end = int(
        (
            int(params["dend_length"] / params["diffusion_length"])
            - (
                int(
                    params["seq_dx"]
                    * params["dispersion_factor"]
                    / params["diffusion_length"]
                )
                * (params["num_inputs"] - 1)
                + 1
            )
        )
        / 2
    )
    for n in range(params["num_inputs"]):
        location = (
            blank_voxels_at_end
            + int(
                (params["seq_dx"] * params["dispersion_factor"])
                / params["diffusion_length"]
            )
            * n
        )
        stim2[params["pre_stim_time"] + n * params["seq_dt"]].append((location, ampl))
        stim2[
            params["pre_stim_time"] + n * params["seq_dt"] + params["stim_width"]
        ].append((location, basal_ampl))
    dispersed_tot = run_stimulus(stim2)

    with open(
        DATA_PATH
        + "sample_dispersed_output_seq_dx_%.6f_seq_dt_%.3f.csv"
        % (params["seq_dx"] * params["dispersion_factor"], params["seq_dt"]),
        "w",
    ) as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(dispersed_tot)


def main():
    """Create a compartment and simulate ca_m activation"""
    global params
    library = moose.Neutral("/library")
    if not os.path.exists(DATA_PATH):
        os.makedirs(DATA_PATH)
    make_passive_soma("cell", params["dend_length"], params["dend_diameter"])
    simulate_ca_cam()
    plt.show()


if __name__ == "__main__":
    main()
