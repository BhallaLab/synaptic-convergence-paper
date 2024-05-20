import numpy as np
import csv
import time
import sys
import os
from numba import njit
import group_params as prm

# np.set_printoptions(threshold=np.inf)


MIN_M = prm.MIN_M
MAX_M = prm.MAX_M

T_POP_POST = 4000

RUN_NUM = int(sys.argv[1])
SEED = int(sys.argv[2])
STIMULUS_SEED = 1947
NUM_TRIALS = 1

DATA_PATH = "./data/"


networks = [
    prm.hippo_chem,
    prm.hippo_cicr,
    # prm.hippo_elec,
    prm.cortex_chem,
    prm.cortex_cicr,
    prm.cortex_elec,
]


def gen_connectivity(t_pop_pre, num_synapses):
    """Generate random connectivity from presynaptic population to postsynaptic population

    Parameters
    ----------
    t_pop_pre : int
        Size of the presynaptic population
    num_synapses : int
        Number of synapses

    Returns
    -------
    pre_conn : np.ndarray
        2D-array of IDs of presynaptic neurons with a size of (T_POP_POST x num_synapses)
    """
    pre_conn = np.random.randint(t_pop_pre, size=(T_POP_POST, num_synapses))
    return pre_conn


def simulate_presynaptic_activity(
    t_pop_pre,
    noise_prob,
    num_ensembles,
    ensemble_size,
    ensemble_participation_prob,
    stim_type,
):
    match stim_type:
        case 0:
            pre_pop_activity = np.random.choice(
                np.arange(-1, 1, dtype=np.int8),
                size=t_pop_pre,
                p=[noise_prob, 1 - noise_prob],
            )
        case 1:
            pre_pop_activity = np.random.choice(
                np.arange(-1, 1, dtype=np.int8),
                size=t_pop_pre,
                p=[noise_prob, 1 - noise_prob],
            )
            pre_pop_activity[0 : num_ensembles * ensemble_size] = np.random.choice(
                np.arange(2, dtype=np.int8),
                size=num_ensembles * ensemble_size,
                p=[1 - ensemble_participation_prob, ensemble_participation_prob],
            )
        case -1:
            pre_pop_activity = np.zeros(t_pop_pre, dtype=np.int8)
            pre_pop_activity[0 : num_ensembles * ensemble_size] = 1
        case _:
            pre_pop_activity = np.zeros(t_pop_pre, dtype=np.int8)
    return pre_pop_activity


def update_presynaptic_activity(
    noise_prob, ensemble_num, ensemble_size, pre_pop_activity, stim_type
):
    match stim_type:
        case 1:
            pre_pop_activity[
                (ensemble_num - 1) * ensemble_size : ensemble_num * ensemble_size
            ] = np.random.choice(
                np.arange(-1, 1, dtype=np.int8),
                size=ensemble_size,
                p=[noise_prob, 1 - noise_prob],
            )
        case -1:
            pre_pop_activity[
                (ensemble_num - 1) * ensemble_size : ensemble_num * ensemble_size
            ] = 0
    return pre_pop_activity


@njit
def generate_input_array(
    pre_connections,
    pre_pop_activity,
    ensemble_size,
    input_array,
    num_ensembles,
    stim_type,
):
    if stim_type != 0:
        for syn, pre_neuron_id in enumerate(pre_connections):
            e_id = (pre_neuron_id // ensemble_size) + 1
            if e_id > num_ensembles:
                input_array[syn] = pre_pop_activity[pre_neuron_id]
            else:
                input_array[syn] = pre_pop_activity[pre_neuron_id] * e_id

        return input_array

    for syn, pre_neuron_id in enumerate(pre_connections):
        input_array[syn] = pre_pop_activity[pre_neuron_id]

    return input_array


@njit
def check_for_groups(
    input_array,
    neuron_id,
    trial_num,
    stimulus_driven_group_splits,
    num_synapses,
    zone_size,
    M,
    stim_type,
):
    for start_ind in range(num_synapses - zone_size + 1):
        if stim_type != 0:
            window = input_array[start_ind : start_ind + zone_size]
            ensemble_input_locs = np.where(window > 0)[0]
            if len(ensemble_input_locs) >= M:
                num_unique_ensembles = len(np.unique(window[ensemble_input_locs]))
                stimulus_driven_group_splits[
                    M - MIN_M, neuron_id, num_unique_ensembles - 1
                ] += 1
    return


def count_groups(t_pop_pre, N, p_e, p_conn, L, Z, R, D):
    stimulus_driven_group_splits = np.zeros(
        (MAX_M - MIN_M + 1, T_POP_POST, MAX_M), dtype=np.uint64
    )

    num_synapses = int(p_conn * t_pop_pre)
    zone_size = int(Z * num_synapses / L)
    noise_prob = 1 - np.exp(-R * D)

    np.random.seed(SEED)
    pre_conn = gen_connectivity(t_pop_pre, num_synapses)

    np.random.seed(STIMULUS_SEED)
    stimulus = np.random.randint(2, size=NUM_TRIALS, dtype=np.int8)
    stimulus[
        0
    ] = (
        -1
    )  # Corresponds to the case where all neurons in ensembles are active, to estimate connectivity based groups

    input_array = np.empty(num_synapses, dtype=np.int32)

    for t_num in range(NUM_TRIALS):
        pre_pop_activity = simulate_presynaptic_activity(
            t_pop_pre, noise_prob, MAX_M, N, p_e, stimulus[t_num]
        )
        for M in range(MAX_M, MIN_M - 1, -1):
            for neuron_id in range(T_POP_POST):
                input_array = generate_input_array(
                    pre_conn[neuron_id],
                    pre_pop_activity,
                    N,
                    input_array,
                    M,
                    stimulus[t_num],
                )
                check_for_groups(
                    input_array,
                    neuron_id,
                    t_num,
                    stimulus_driven_group_splits,
                    num_synapses,
                    zone_size,
                    M,
                    stimulus[t_num],
                )
            pre_pop_activity = update_presynaptic_activity(
                noise_prob, M, N, pre_pop_activity, stimulus[t_num]
            )

    return (
        stimulus,
        stimulus_driven_group_splits,
    )


def process_network(network):
    (
        stimulus,
        stimulus_driven_pos_group_splits,
    ) = count_groups(
        t_pop_pre=network["population_size"],
        N=network["ensemble_size"],
        p_e=network["ensemble_participation_prob"],
        p_conn=network["connection_prob"],
        L=network["dendrite_length"],
        Z=network["zone_length"],
        R=network["background_rate"],
        D=network["input_duration"],
    )

    for M in range(MIN_M, MAX_M + 1):
        filename = DATA_PATH + "stimulus-driven-group-splits-%s-run-%03d-M-%d.csv" % (
            network["name"],
            RUN_NUM,
            M,
        )
        with open(filename, "w") as f:
            writer = csv.writer(f, delimiter=",")
            for row in stimulus_driven_pos_group_splits[M - MIN_M]:
                writer.writerow(row[:M])

    return


def main():
    t0 = time.time()
    if not os.path.exists(DATA_PATH):
        os.makedirs(DATA_PATH)
    for net in networks:
        print(f"Run num = {RUN_NUM}, seed = {SEED}, net = {net['name']}")
        process_network(net)
    print("runtime=", time.time() - t0)


if __name__ == "__main__":
    main()
