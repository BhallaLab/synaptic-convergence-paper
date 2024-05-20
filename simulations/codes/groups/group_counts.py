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
NUM_TRIALS = 100

DATA_PATH = "./data/"


networks = [
    prm.hippo_chem,
    prm.hippo_cicr,
    prm.hippo_elec,
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
        # print("Stim type is 1")
    # print(stim_type)
    # print(pre_pop_activity[100:])
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
    # print(pre_pop_activity[100:])
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
    # print(f"{pre_connections[:100]=}")
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
    fully_mixed_groups,
    stimulus_driven_groups,
    noise_groups,
    all_groups,
    num_synapses,
    zone_size,
    M,
    stim_type,
):
    for start_ind in range(num_synapses - zone_size + 1):
        if stim_type != 0:
            fully_mixed_group_found = True

            for m in range(1, M + 1):
                if m not in input_array[start_ind : start_ind + zone_size]:
                    fully_mixed_group_found = False
                    break
            if fully_mixed_group_found:
                fully_mixed_groups[trial_num, neuron_id, M - MIN_M] += 1

            if (
                len(np.where(input_array[start_ind : start_ind + zone_size] > 0)[0])
                >= M
            ):
                stimulus_driven_groups[trial_num, neuron_id, M - MIN_M] += 1

            if (
                len(np.where(input_array[start_ind : start_ind + zone_size] != 0)[0])
                >= M
            ):
                all_groups[trial_num, neuron_id, M - MIN_M] += 1

        if stim_type != -1:
            if (
                len(np.where(input_array[start_ind : start_ind + zone_size] < 0)[0])
                >= M
            ):
                noise_groups[trial_num, neuron_id, M - MIN_M] += 1

    return


def count_groups(t_pop_pre, N, p_e, p_conn, L, Z, R, D):
    fully_mixed_groups = np.zeros(
        (NUM_TRIALS, T_POP_POST, MAX_M - MIN_M + 1), dtype=np.uint64
    )
    stimulus_driven_groups = np.zeros(
        (NUM_TRIALS, T_POP_POST, MAX_M - MIN_M + 1), dtype=np.uint64
    )
    noise_groups = np.zeros(
        (NUM_TRIALS, T_POP_POST, MAX_M - MIN_M + 1), dtype=np.uint64
    )
    all_groups = np.zeros((NUM_TRIALS, T_POP_POST, MAX_M - MIN_M + 1), dtype=np.uint64)

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

    # print(f"{pre_conn=}\n{pre_conn.shape}")
    input_array = np.empty(num_synapses, dtype=np.int32)

    for t_num in range(NUM_TRIALS):
        # print(f"{t_num=}, {stimulus[t_num]=}")
        pre_pop_activity = simulate_presynaptic_activity(
            t_pop_pre, noise_prob, MAX_M, N, p_e, stimulus[t_num]
        )
        # print(f"{pre_pop_activity[:1000]}\n{pre_pop_activity.shape}")
        for M in range(MAX_M, MIN_M - 1, -1):
            # print(f"{M=}")
            for neuron_id in range(T_POP_POST):
                # print(f"{neuron_id=}")
                input_array = generate_input_array(
                    pre_conn[neuron_id],
                    pre_pop_activity,
                    N,
                    input_array,
                    M,
                    stimulus[t_num],
                )
                # print(f"{input_array[:100]}\n {input_array.shape}")
                check_for_groups(
                    input_array,
                    neuron_id,
                    t_num,
                    fully_mixed_groups,
                    stimulus_driven_groups,
                    noise_groups,
                    all_groups,
                    num_synapses,
                    zone_size,
                    M,
                    stimulus[t_num],
                )
            pre_pop_activity = update_presynaptic_activity(
                noise_prob, M, N, pre_pop_activity, stimulus[t_num]
            )
            # print(f"{pre_pop_activity[:1000]}\n{pre_pop_activity.shape}")

    return (
        stimulus,
        fully_mixed_groups,
        stimulus_driven_groups,
        noise_groups,
        all_groups,
    )


def process_network(network):
    (
        stimulus,
        fully_mixed_pos_groups,
        stimulus_driven_pos_groups,
        noise_pos_groups,
        all_pos_groups,
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

    filename = DATA_PATH + "stimulus-%s-run-%03d.csv" % (network["name"], RUN_NUM)
    with open(filename, "w") as f:
        writer = csv.writer(f, delimiter=",")
        writer.writerow(stimulus)

    for t_num in range(NUM_TRIALS):
        filename = DATA_PATH + "noise-groups-%s-run-%03d-trial-%03d.csv" % (
            network["name"],
            RUN_NUM,
            t_num + 1,
        )
        with open(filename, "w") as f:
            writer = csv.writer(f, delimiter=",")
            for row in noise_pos_groups[t_num]:
                writer.writerow(row)

        if stimulus[t_num]:
            filename = DATA_PATH + "fully-mixed-groups-%s-run-%03d-trial-%03d.csv" % (
                network["name"],
                RUN_NUM,
                t_num + 1,
            )
            with open(filename, "w") as f:
                writer = csv.writer(f, delimiter=",")
                for row in fully_mixed_pos_groups[t_num]:
                    writer.writerow(row)

            filename = (
                DATA_PATH
                + "stimulus-driven-groups-%s-run-%03d-trial-%03d.csv"
                % (
                    network["name"],
                    RUN_NUM,
                    t_num + 1,
                )
            )
            with open(filename, "w") as f:
                writer = csv.writer(f, delimiter=",")
                for row in stimulus_driven_pos_groups[t_num]:
                    writer.writerow(row)

            filename = DATA_PATH + "all-groups-%s-run-%03d-trial-%03d.csv" % (
                network["name"],
                RUN_NUM,
                t_num + 1,
            )
            with open(filename, "w") as f:
                writer = csv.writer(f, delimiter=",")
                for row in all_pos_groups[t_num]:
                    writer.writerow(row)

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
