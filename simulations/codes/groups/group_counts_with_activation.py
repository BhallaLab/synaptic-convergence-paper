import numpy as np
import csv
import time
import sys
import os
from numba import njit
import group_params as prm


MIN_M = prm.MIN_M
MAX_M = 4
NUM_POWERS = 5

T_POP_POST = 4000

RUN_NUM = int(sys.argv[1])
SEED = int(sys.argv[2])
STIMULUS_SEED = 1947
NUM_TRIALS = 200

DATA_PATH = "./activation_data/"


networks = [
    prm.hippo_elec,
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
def calc_power_activation(
    input_array, neuron_id, trial_num, activation, num_synapses, zone_size, M
):
    zonal_input = [
        sum(input_array[i : i + zone_size] != 0) for i in range(num_synapses)
    ]
    filtered_actn = [1] * num_synapses
    powers = range(1, NUM_POWERS + 1)
    for power in powers:
        filtered_actn = [filtered_actn[i] * zonal_input[i] for i in range(num_synapses)]
        activation[trial_num, power - 1, neuron_id, M - MIN_M] = (
            sum(filtered_actn) / zone_size
        )

    return activation, powers


@njit
def calc_weibull_exponential_activation(
    input_array, neuron_id, trial_num, activation, num_synapses, zone_size, M
):
    zonal_input = [
        sum(input_array[i : i + zone_size] != 0) for i in range(num_synapses)
    ]
    zonal_input = np.array(zonal_input)

    complete_group_actn = (
        0.9  # Activation achieved when all M inputs are present. Must be between 0 to 1
    )
    # First we will calulcate the linear summation. Note that this is different from the weibull based function. We need it for a comparision, hence doing it in the same place.
    # Complete group activation = 0.9. So the function needs to produce a value of 0.9 when there are M inputs. Ultimately, the function is capped at 1
    # For the linear summation case, y=kx, If y=0.9 at x=M, k = 0.9/M
    c = complete_group_actn / M
    filtered_actn = [c * zonal_input[i] for i in range(num_synapses)]
    activation[trial_num, 0, neuron_id, M - MIN_M] = sum(filtered_actn) / zone_size

    # Now we begin the Weibull based calculations
    ratio = 2
    lambda_val = M
    base = 1 / (1 - complete_group_actn)
    exponents = np.array([ratio**i for i in range(NUM_POWERS)])
    filtered_actn = [1] * num_synapses
    for p, k in enumerate(exponents[1:]):
        filtered_actn = [1 - base ** (-((x / lambda_val) ** k)) for x in zonal_input]
        activation[trial_num, p + 1, neuron_id, M - MIN_M] = (
            sum(filtered_actn) / zone_size
        )

    return activation, exponents


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


def count_groups(t_pop_pre, N, p_e, p_conn, L, Z, R, D, nonlinearity):
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

    activation = np.zeros(
        (NUM_TRIALS, NUM_POWERS + 1, T_POP_POST, MAX_M - MIN_M + 1),
        dtype=np.double,
    )
    match nonlinearity:
        case "power":
            calc_activation = calc_power_activation
        case "weibull_exponent":
            calc_activation = calc_weibull_exponential_activation

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
                activation, nonlinearity_parameters = calc_activation(
                    input_array,
                    neuron_id,
                    t_num,
                    activation,
                    num_synapses,
                    zone_size,
                    M,
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
        activation,
        nonlinearity_parameters,
    )


def process_network(network):
    nonlinearity = "weibull_exponent"
    (
        stimulus,
        fully_mixed_pos_groups,
        stimulus_driven_pos_groups,
        noise_pos_groups,
        all_pos_groups,
        activation_matrix,
        nonlinearity_parameters,
    ) = count_groups(
        t_pop_pre=network["population_size"],
        N=network["ensemble_size"],
        p_e=network["ensemble_participation_prob"],
        p_conn=network["connection_prob"],
        L=network["dendrite_length"],
        Z=network["zone_length"],
        R=network["background_rate"],
        D=network["input_duration"],
        nonlinearity=nonlinearity,
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

            for p, param in enumerate(nonlinearity_parameters):
                filename = (
                    DATA_PATH
                    + "activation-of-groups-%s-run-%03d-trial-%03d-%s-%d.csv"
                    % (
                        network["name"],
                        RUN_NUM,
                        t_num + 1,
                        nonlinearity,
                        param,
                    )
                )
                with open(filename, "w") as f:
                    writer = csv.writer(f, delimiter=",")
                    for row in activation_matrix[t_num, p]:
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
