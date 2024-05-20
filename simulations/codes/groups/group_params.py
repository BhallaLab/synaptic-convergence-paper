MIN_M = 3
MAX_M = 9
MIN_POWER = 1
MAX_POWER = 5


hippo_chem = {
    "name": "hippo-chem",
    "population_size": 400000,
    "connection_prob": 0.05,
    "ensemble_size": 100,
    "ensemble_participation_prob": 0.8,
    "zone_length": 10e-6,
    "dendrite_length": 10e-3,
    "background_rate": 1e-2,
    "input_duration": 2,
    "spine_spacing": 5e-7,
}

hippo_cicr = {
    "name": "hippo-CICR",
    "population_size": 400000,
    "connection_prob": 0.05,
    "ensemble_size": 100,
    "ensemble_participation_prob": 0.8,
    "zone_length": 10e-6,
    "dendrite_length": 10e-3,
    "background_rate": 1e-2,
    "input_duration": 2e-1,
    "spine_spacing": 5e-7,
}

hippo_elec = {
    "name": "hippo-elec",
    "population_size": 400000,
    "connection_prob": 0.05,
    "ensemble_size": 100,
    "ensemble_participation_prob": 0.8,
    "zone_length": 50e-6,
    "dendrite_length": 10e-3,
    "background_rate": 1e-1,
    "input_duration": 4e-3,
    "spine_spacing": 5e-7,
}

cortex_chem = {
    "name": "cortex-chem",
    "population_size": 100000,
    "connection_prob": 0.2,
    "ensemble_size": 100,
    "ensemble_participation_prob": 0.8,
    "zone_length": 10e-6,
    "dendrite_length": 10e-3,
    "background_rate": 1e-1,
    "input_duration": 2,
    "spine_spacing": 5e-7,
}

cortex_cicr = {
    "name": "cortex-CICR",
    "population_size": 100000,
    "connection_prob": 0.2,
    "ensemble_size": 100,
    "ensemble_participation_prob": 0.8,
    "zone_length": 10e-6,
    "dendrite_length": 10e-3,
    "background_rate": 1e-1,
    "input_duration": 2e-1,
    "spine_spacing": 5e-7,
}

cortex_elec = {
    "name": "cortex-elec",
    "population_size": 100000,
    "connection_prob": 0.2,
    "ensemble_size": 100,
    "ensemble_participation_prob": 0.8,
    "zone_length": 50e-6,
    "dendrite_length": 10e-3,
    "background_rate": 1.0,
    "input_duration": 4e-3,
    "spine_spacing": 5e-7,
}
