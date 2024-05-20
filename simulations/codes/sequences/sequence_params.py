MIN_M = 3
MAX_M = 9

T_POP_POST = 10

hippo_chem = {
    "name": "hippo-chem",
    "population_size": 400000,
    "connection_prob": 0.05,
    "ensemble_size": 1000,
    "ensemble_participation_prob": 0.8,
    "input_spacing": 2e-6,
    "dendrite_length": 10e-3,
    "background_rate": 1e-2,
    "sequence_time_step": 2.0,
    "flank_length": 0.0,
    "spine_spacing": 5e-7,
    "delta": 1.5e-6,
    "decay_factor": 0.8,
    "nonlinear_factor": 1.2,
    "sigmoid_exponent": 0.1,
    "sigmoid_max_value": 80,
}

hippo_cicr = {
    "name": "hippo-CICR",
    "population_size": 400000,
    "connection_prob": 0.05,
    "ensemble_size": 1000,
    "ensemble_participation_prob": 0.8,
    "input_spacing": 2e-6,
    "dendrite_length": 10e-3,
    "background_rate": 1e-2,
    "sequence_time_step": 2e-1,
    "flank_length": 0.0,
    "spine_spacing": 5e-7,
    "delta": 1.5e-6,
    "decay_factor": 0.8,
    "nonlinear_factor": 1.2,
    "sigmoid_exponent": 0.1,
    "sigmoid_max_value": 80,
}

hippo_elec = {
    "name": "hippo-elec",
    "population_size": 400000,
    "connection_prob": 0.05,
    "ensemble_size": 1000,
    "ensemble_participation_prob": 0.8,
    "input_spacing": 10e-6,
    "dendrite_length": 10e-3,
    "background_rate": 1e-1,
    "sequence_time_step": 4e-3,
    "flank_length": 100e-6,
    "spine_spacing": 5e-7,
    "delta": 5e-6,
    "decay_factor": 0.85,  # Assuming a tau decay of 25ms, and time of 4ms
    "nonlinear_factor": 0.4,  # To obtain a selectivity of ~10%
}

cortex_chem = {
    "name": "cortex-chem",
    "population_size": 100000,
    "connection_prob": 0.2,
    "ensemble_size": 1000,
    "ensemble_participation_prob": 0.8,
    "input_spacing": 2e-6,
    "dendrite_length": 10e-3,
    "background_rate": 1e-1,
    "sequence_time_step": 2.0,
    "flank_length": 0.0,
    "spine_spacing": 5e-7,
    "delta": 1.5e-6,
    "decay_factor": 0.8,
    "nonlinear_factor": 1.2,
    "sigmoid_exponent": 0.1,
    "sigmoid_max_value": 80,
}

cortex_cicr = {
    "name": "cortex-CICR",
    "population_size": 100000,
    "connection_prob": 0.2,
    "ensemble_size": 1000,
    "ensemble_participation_prob": 0.8,
    "input_spacing": 2e-6,
    "dendrite_length": 10e-3,
    "background_rate": 1e-1,
    "sequence_time_step": 2e-1,
    "flank_length": 0.0,
    "spine_spacing": 5e-7,
    "delta": 1.5e-6,
    "decay_factor": 0.8,
    "nonlinear_factor": 1.2,
    "sigmoid_exponent": 0.1,
    "sigmoid_max_value": 80,
}

cortex_elec = {
    "name": "cortex-elec",
    "population_size": 100000,
    "connection_prob": 0.2,
    "ensemble_size": 1000,
    "ensemble_participation_prob": 0.8,
    "input_spacing": 10e-6,
    "dendrite_length": 10e-3,
    "background_rate": 1.0,
    "sequence_time_step": 4e-3,
    "flank_length": 100e-6,
    "spine_spacing": 5e-7,
    "delta": 5e-6,
    "decay_factor": 0.85,  # Assuming a tau decay of 25ms, and time of 4ms
    "nonlinear_factor": 0.4,  # To obtain a selectivity of ~10%
}
