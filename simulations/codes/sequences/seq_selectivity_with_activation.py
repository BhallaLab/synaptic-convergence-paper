from numba.parfors.parfor import parfor_add_offset_to_labels
import numpy as np
import csv
import time
import sys
import sequence_params as prm
import os
from numba import njit
from scipy import sparse
from itertools import permutations


T_POP_POST = 10000
RUN_NUM = int(sys.argv[1])
SEED = int(sys.argv[2])
STIMULUS_SEED = 47
NUM_TRIALS = 100
PATTERN_LENGTH = 4
MAX_NUM_PATTERNS = 24

DATA_PATH = "./activation_data/"

networks = [prm.hippo_cicr]


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
    stim_order,
    pre_stim_period,
    post_stim_period,
    ensemble_size,
    ensemble_participation_prob,
    stim_type,
):
    """Simulate presynaptic activity for the duration of a sequence

    Parameters
    ----------
    t_pop_pre : int
        Size of the presynaptic population
    noise_prob : float
        Probability of background activity
    stim_order : array like
        Order of stimulation of the ensembles. [3,1,2,4,5] implies ensemble #3 is activated first, followed by ensemble #1 and so on.
    pre_stim_period : int
        Number of extra time steps to stimulate background activity before the start of the sequence
    post_stim_period : int
        Number of extra time steps to stimulate background activity while waiting for the sequence activations to decay
    ensemble_size : int
        Size of the ensemble in terms of the number of neurons
    ensemble_participation_prob : float
        Probability of an ensemble neuron being active
    stim_type : int
        Nature of the stimulus
        -1: All ensemble neurons are active due to the stimulus
        0: No ensemble activity, only background activity
        1: Ensembles with a fraction of neurons active, along with background activity in the non-ensemblar neurons

    Returns
    -------
    pre_pop_activity : np.ndarray
        2D array of the size of (t_pop_pre x len(stim_order)), representing presynaptic activity.
        -1: Active due to background activity
        0: Inactive
        1: Active due to participation in the ensemble
    """
    total_time = pre_stim_period + len(stim_order) + post_stim_period
    match stim_type:
        case 0:
            pre_pop_activity = np.random.choice(
                np.arange(-1, 1, dtype=np.int8),
                size=(t_pop_pre, total_time),
                p=[noise_prob, 1 - noise_prob],
            )
        case 1:
            pre_pop_activity = np.random.choice(
                np.arange(-1, 1, dtype=np.int8),
                size=(t_pop_pre, total_time),
                p=[noise_prob, 1 - noise_prob],
            )
            for i, ensemble_id in enumerate(stim_order):
                pre_pop_activity[
                    (ensemble_id - 1) * ensemble_size : ensemble_id * ensemble_size,
                    i + pre_stim_period,
                ] = np.random.choice(
                    [1, 0],
                    size=ensemble_size,
                    p=[ensemble_participation_prob, 1 - ensemble_participation_prob],
                )
        case -1:
            pre_pop_activity = np.zeros((t_pop_pre, total_time), dtype=np.int8)
            for i, ensemble_id in enumerate(stim_order):
                pre_pop_activity[
                    (ensemble_id - 1) * ensemble_size : ensemble_id * ensemble_size,
                    i + pre_stim_period,
                ] = 1
        case _:
            pre_pop_activity = np.zeros((t_pop_pre, total_time), dtype=np.int8)
    return pre_pop_activity


@njit
def update_ensemble_activity(
    pre_pop_activity,
    stim_order,
    prev_stim_order,
    ensemble_size,
    pre_stim_period,
):
    """Update the activity of ensembles based on order of stimulus

    Parameters
    ----------
    pre_pop_activity : np.ndarray
        2D array of the size of (t_pop_pre x len(stim_order)), representing presynaptic activity.
    stim_order : array like
        Order of stimulation of the ensembles. [3,1,2,4,5] implies ensemble #3 is activated first, followed by ensemble #1 and so on.

    prev_stim_order : array like
        Order of stimulation of the ensembles in the previous pattern
    ensemble_size : int
        Size of the ensemble in terms of the number of neurons
    pre_stim_period : int
        Number of time steps before the onset of ensemble activity

    Returns
    -------
    pre_pop_activity : np.ndarray
        2D array of the size of (t_pop_pre x len(stim_order)), representing the updated presynaptic activity.
        -1: Active due to background activity
        0: Inactive
        1: Active due to participation in the ensemble
    """
    for i, ensemble_id in enumerate(stim_order):
        prev_stim_time = (
            np.where(prev_stim_order == ensemble_id)[0][0] + pre_stim_period
        )  # When was this ensemble active in the previous stim order
        temp_activity = np.copy(
            pre_pop_activity[
                (ensemble_id - 1) * ensemble_size : ensemble_id * ensemble_size,
                pre_stim_period + i,
            ]
        )
        pre_pop_activity[
            (ensemble_id - 1) * ensemble_size : ensemble_id * ensemble_size,
            pre_stim_period + i,
        ] = np.copy(
            pre_pop_activity[
                (ensemble_id - 1) * ensemble_size : ensemble_id * ensemble_size,
                prev_stim_time,
            ]
        )
        pre_pop_activity[
            (ensemble_id - 1) * ensemble_size : ensemble_id * ensemble_size,
            prev_stim_time,
        ] = np.copy(temp_activity)
    return pre_pop_activity


@njit
def gen_dist_matrix(
    seq_type,
    prev_activity,
    curr_activity,
    pre_IDs,
    num_synapses,
    spacing_size,
    input_width,
):
    """Generate distance based activity relation matrix from activity at the previous and current time steps. If
    synapse i and j received input at t-D and t respectively, and spacing_size < j-i < spacing_size+input_width, the
    i,j th entry in the matrix is set to 1

    Parameters
    ----------
    seq_type : str
        Type of sequence: poss/noise/all
    prev_activity : np.array
        Activity vector of the pre-synaptic population at t-D
    curr_activity : np.array
        Activity vectorof the presynaptic population at t
    pre_IDs : np.array
        Vector of presynaptic neuron IDs connecting to the synapses
    num_synapses : int
        Number of synapses on a postsynaptic neuron
    spacing_size : int
        Minimum number of synapses that constitute the gap between two sequential inputs; derived from spacing window S
    input_width : int
        Number of synapses that fall within the spacing considered for sequences to occur; derived from delta

    Returns
    -------
    (data, (row_ind, col_ind)) : tuple
        Tuple of the data entries, row indices and column indices needed to construct a sparse matrix
    """
    data = []
    row_ind = []
    col_ind = []
    match seq_type:
        case "poss":
            for i in range(num_synapses):
                for j in range(
                    min(i + spacing_size, num_synapses),
                    min(i + spacing_size + input_width, num_synapses),
                ):
                    if (
                        prev_activity[pre_IDs[i]] == 1
                        and curr_activity[pre_IDs[j]] == 1
                    ):
                        row_ind.append(i)
                        col_ind.append(j)
                        data.append(1)
        case "noise":
            for i in range(num_synapses):
                for j in range(
                    min(i + spacing_size, num_synapses),
                    min(i + spacing_size + input_width, num_synapses),
                ):
                    if (
                        prev_activity[pre_IDs[i]] == -1
                        and curr_activity[pre_IDs[j]] == -1
                    ):
                        row_ind.append(i)
                        col_ind.append(j)
                        data.append(1)
        case "all":
            for i in range(num_synapses):
                for j in range(
                    min(i + spacing_size, num_synapses),
                    min(i + spacing_size + input_width, num_synapses),
                ):
                    if (
                        prev_activity[pre_IDs[i]] != 0
                        and curr_activity[pre_IDs[j]] != 0
                    ):
                        row_ind.append(i)
                        col_ind.append(j)
                        data.append(1)
    return (data, (row_ind, col_ind))


def mul_mat(A, B):
    """Multiply two matrices

    Parameters
    ----------
    A : sparse.csr_array
        Sparse matrix 1
    B : sparse.csr_array
        Sparse matrix 2

    Returns
    -------
    C : sparse.csr_array
        Product of the two matrices
    """
    C = A.dot(B)
    return C


@njit
def generate_dendritic_input_matrix(
    pre_connections,
    pre_pop_activity,
    total_time,
):
    """Generate dendritic input matrix

    Parameters
    ----------
    pre_connections : np.ndarray
        1D-array of IDs of presynaptic neurons with a size of num_synapses
    pre_pop_activity : np.ndarray
        2D array of the size of (t_pop_pre x len(stim_order)), representing presynaptic activity.
        -1: Active due to background activity
        0: Inactive
        1: Active due to participation in the ensemble
    total_time : int
        Total duration in time steps

    Returns
    -------
    input_matrix : np.ndarray
        (num_synapses x seq_len) array containing 1s where inputs were active at tha time point and 0s for inactive inputs
    """
    input_matrix = np.empty((total_time, len(pre_connections)), dtype=np.int32)
    for m in range(total_time):
        for syn, pre_neuron_id in enumerate(pre_connections):
            input_matrix[m, syn] = abs(pre_pop_activity[pre_neuron_id, m])

    return input_matrix


@njit
def semi_sigmoid(x, max_val, exponent):
    """Return output from a y-shifted sigmoid, whose midpoint in the y-axis passes through (0,0)

    Parameters
    ----------
    x : float
        x value
    max_val : float
        Saturation value of the sigmoid funtion as x-> infinity
    exponent : float
        Exponent value that controls the steepness of the sigmoid

    Returns
    -------
    output : float
        Output from the y-shifted sigmoid function
    """
    output = max_val * 2 * ((1 / (1 + np.exp(-exponent * x))) - 0.5)
    return output


@njit
def get_activation_peaks_and_auc(
    dendritic_input_matrix,
    spacing_size,
    input_width,
    gamma,
    eta,
    max_val,
    exponent,
):
    """Get activation peaks and area under the curve

    Parameters
    ----------
    dendritic_input_matrix : np.ndarray
        ()
    spacing_size : int
        Spacing size based on S, in number of synapses
    input_width : int
        Input width (based on delta) in number of synapses
    gamma : float
        Decay factor
    eta : float
        Nonlinear factor
    max_val : float
        Saturation value of the sigmoid funtion as x-> infinity
    exponent : float
        Exponent value that controls the steepness of the sigmoid

    Returns
    -------
    peak : float
        Peak value activation summed over the entire dendrite
    auc : float
        Area under the curve of activation summed over the entire dendrite over the entire time window
    """
    activation = np.zeros_like(
        dendritic_input_matrix,
        dtype=np.float64,
    )
    activation[0, :] = dendritic_input_matrix[0, :]

    for t in range(1, dendritic_input_matrix.shape[0]):
        for syn in range(dendritic_input_matrix.shape[1]):
            # Existing activation decays by gamma in each time step
            activation[t, syn] = gamma * activation[t - 1, syn]
            if dendritic_input_matrix[t, syn]:
                # If an input arrives, it has a baseline contribution and a sequential contribution
                activation[t, syn] += dendritic_input_matrix[t, syn] + semi_sigmoid(
                    sum(
                        [
                            (
                                1
                                + gamma
                                * activation[t - 1, j]
                                * dendritic_input_matrix[t, syn]
                            )
                            ** eta
                            - 1
                            for j in range(
                                max(syn - spacing_size - input_width + 1, 0),
                                max(0, syn - spacing_size + 1),
                            )
                        ]
                    ),
                    max_val=max_val,
                    exponent=exponent,
                )

    actn_time_profile = np.sum(activation, axis=1)
    peak = np.max(actn_time_profile)
    auc = np.sum(actn_time_profile)
    return peak, auc


def count_sequences(
    t_pop_pre, N, p_conn, p_e, L, S, delta, R, D, gamma, eta, max_val, exponent
):
    """Count different kinds of sequences occurring on a neuron's dendrite, obtain the peaks and area under the curve

    Parameters
    ----------
    t_pop_pre : int
        Size of the presynaptic population
    N : int
        Ensemble size
    p_conn : float
        Probability of connectivity between the presynaptic and postsynaptic neurons
    p_e : float
        Probability of an ensemble neuron being active, i.e., the probability of it's participation in the ensemble
    L : float
        Length of the dendrite
    S : float
        Minimum spacing between two inputs
    delta : float
        Available zone length for an input to be a part of a sequence
    R : float
        Rate of background activity in the network
    D : float
        Sequence time step
    gamma : float
        Activation decay factor
    eta : float
        Nonlinear factor
    max_val : float
        Saturation value of the sigmoid funtion as x-> infinity
    exponent : float
        Exponent value that controls the steepness of the sigmoid

    Returns
    -------
    stimulus : np.ndarray
        1D-array of stimulus type of size number of trials containing values from [-1, 0, 1].
        -1: All ensemble neurons are active due to the stimulus
        0: No ensemble activity, only background activity
        1: Ensembles with a fraction of neurons active, along with background activity in the non-ensemblar neurons
    poss_sequences : np.ndarray
        Array of size (num_trials x postsyn_population_size x num_seq_lengths) containing the number of poss sequences on a neuron
    dend_activation_peak : np.ndarray
        Array of size (num_trials x postsyn_population_size x num_seq_lengths) containing the peak activation of a neuron
    dend_activation_auc : np.ndarray
        Array of size (num_trials x postsyn_population_size x num_seq_lengths) containing the total area under the curve activation of a neuron
    """
    num_synapses = int(p_conn * t_pop_pre)
    spacing_size = int(S * num_synapses / L)
    input_width = int(delta * num_synapses / L)
    noise_prob = 1 - np.exp(-R * D)
    post_stim_period = round(np.log(0.1) / np.log(gamma))
    pre_stim_period = 0

    perm = list(np.array(pat) for pat in permutations(range(1, PATTERN_LENGTH + 1)))
    if len(perm) > MAX_NUM_PATTERNS:
        input_patterns = perm[0 :: int(len(perm) / MAX_NUM_PATTERNS)]
    else:
        input_patterns = perm
    input_patterns.insert(0, input_patterns[0])

    poss_sequences = np.zeros(
        (NUM_TRIALS, T_POP_POST, len(input_patterns) - 1), dtype=np.uint64
    )
    all_sequences = np.zeros(
        (NUM_TRIALS, T_POP_POST, len(input_patterns) - 1), dtype=np.uint64
    )
    dend_activation_peak = np.zeros(
        (NUM_TRIALS, T_POP_POST, len(input_patterns) - 1), dtype=np.float64
    )
    dend_activation_auc = np.zeros(
        (NUM_TRIALS, T_POP_POST, len(input_patterns) - 1), dtype=np.float64
    )

    np.random.seed(SEED)
    pre_conn = gen_connectivity(t_pop_pre, num_synapses)

    np.random.seed(STIMULUS_SEED)
    stimulus = np.random.randint(2, size=NUM_TRIALS, dtype=np.int8)
    stimulus[
        0
    ] = (
        -1
    )  # Corresponds to the case where all neurons in ensembles are active, to estimate connectivity based groups

    for t_num in range(NUM_TRIALS):
        activity = simulate_presynaptic_activity(
            t_pop_pre,
            noise_prob,
            input_patterns[0],
            pre_stim_period,
            post_stim_period,
            N,
            p_e,
            stimulus[t_num],
        )
        if stimulus[t_num]:
            patterns = input_patterns
        else:
            patterns = [
                np.zeros(PATTERN_LENGTH)
            ] * 2  # Input patterns do not have any significance in the case of background trials
        prev_pattern = patterns[0]

        for sp, stim_pattern in enumerate(patterns[1:]):
            if stimulus[t_num]:
                activity = update_ensemble_activity(
                    activity, stim_pattern, prev_pattern, N, pre_stim_period
                )
            for neuron_id in range(T_POP_POST):
                curr_dist_matrix_poss_seq = sparse.csr_array(
                    gen_dist_matrix(
                        "poss",
                        activity[:, 0 + pre_stim_period],
                        activity[:, 1 + pre_stim_period],
                        pre_conn[neuron_id],
                        num_synapses,
                        spacing_size,
                        input_width,
                    ),
                    shape=(num_synapses, num_synapses),
                )
                cum_paths_poss_seq = curr_dist_matrix_poss_seq

                curr_dist_matrix_all_seq = sparse.csr_array(
                    gen_dist_matrix(
                        "all",
                        activity[:, 0 + pre_stim_period],
                        activity[:, 1 + pre_stim_period],
                        pre_conn[neuron_id],
                        num_synapses,
                        spacing_size,
                        input_width,
                    ),
                    shape=(num_synapses, num_synapses),
                )
                cum_paths_all_seq = curr_dist_matrix_all_seq

                for M in range(3, PATTERN_LENGTH + 1):
                    if cum_paths_poss_seq.count_nonzero() > 0:
                        next_dist_matrix_poss_seq = sparse.csr_array(
                            gen_dist_matrix(
                                "poss",
                                activity[:, M - 2 + pre_stim_period],
                                activity[:, M - 1 + pre_stim_period],
                                pre_conn[neuron_id],
                                num_synapses,
                                spacing_size,
                                input_width,
                            ),
                            shape=(num_synapses, num_synapses),
                        )
                        cum_paths_poss_seq = mul_mat(
                            cum_paths_poss_seq, next_dist_matrix_poss_seq
                        )

                    if cum_paths_all_seq.count_nonzero() > 0:
                        next_dist_matrix_all_seq = sparse.csr_array(
                            gen_dist_matrix(
                                "all",
                                activity[:, M - 2 + pre_stim_period],
                                activity[:, M - 1 + pre_stim_period],
                                pre_conn[neuron_id],
                                num_synapses,
                                spacing_size,
                                input_width,
                            ),
                            shape=(num_synapses, num_synapses),
                        )
                        cum_paths_all_seq = mul_mat(
                            cum_paths_all_seq, next_dist_matrix_all_seq
                        )

                poss_sequences[t_num, neuron_id, sp] = np.sum(cum_paths_poss_seq)
                all_sequences[t_num, neuron_id, sp] = np.sum(cum_paths_all_seq)

                dendritic_input_matrix = generate_dendritic_input_matrix(
                    pre_conn[neuron_id],
                    activity,
                    pre_stim_period + PATTERN_LENGTH + post_stim_period,
                )
                peak, auc = get_activation_peaks_and_auc(
                    dendritic_input_matrix,
                    spacing_size,
                    input_width,
                    gamma,
                    eta,
                    max_val,
                    exponent,
                )
                dend_activation_peak[t_num, neuron_id, sp] = peak
                dend_activation_auc[t_num, neuron_id, sp] = auc
            prev_pattern = stim_pattern

    return (
        stimulus,
        input_patterns[1:],
        poss_sequences,
        all_sequences,
        dend_activation_peak,
        dend_activation_auc,
    )


def process_network(network):
    """Simulate presynaptic activity and count sequences for the network

    Parameters
    ----------
    network : dict
        Dictionary containing the parameters of the network
    """
    (
        stimulus,
        patterns,
        poss_pos_sequences,
        all_pos_sequences,
        activation_peak,
        activation_auc,
    ) = count_sequences(
        t_pop_pre=network["population_size"],
        N=network["ensemble_size"],
        p_conn=network["connection_prob"],
        p_e=network["ensemble_participation_prob"],
        L=network["dendrite_length"],
        S=network["input_spacing"],
        delta=network["delta"],
        R=network["background_rate"],
        D=network["sequence_time_step"],
        gamma=network["decay_factor"],
        eta=network["nonlinear_factor"],
        max_val=network["sigmoid_max_value"],
        exponent=network["sigmoid_exponent"],
    )
    filename = DATA_PATH + "patterns-stimulus-%s-run-%03d-M-%d.csv" % (
        network["name"],
        RUN_NUM,
        PATTERN_LENGTH,
    )
    with open(filename, "w") as f:
        writer = csv.writer(f, delimiter=",")
        writer.writerow(stimulus)

    filename = DATA_PATH + "patterns-%s-run-%03d-M-%d.csv" % (
        network["name"],
        RUN_NUM,
        PATTERN_LENGTH,
    )
    with open(filename, "w") as f:
        writer = csv.writer(f, delimiter=",")
        for pattern in patterns:
            writer.writerow(pattern)

    for t_num in range(NUM_TRIALS):
        if stimulus[t_num]:
            pattern_limit = MAX_NUM_PATTERNS
            filename = (
                DATA_PATH
                + "patterns-poss-sequences-%s-run-%03d-trial-%03d-M-%d.csv"
                % (
                    network["name"],
                    RUN_NUM,
                    t_num + 1,
                    PATTERN_LENGTH,
                )
            )
            with open(filename, "w") as f:
                writer = csv.writer(f, delimiter=",")
                for row in poss_pos_sequences[t_num]:
                    writer.writerow(row[:pattern_limit])

        else:
            pattern_limit = 1

        filename = (
            DATA_PATH
            + "patterns-all-sequences-%s-run-%03d-trial-%03d-M-%d.csv"
            % (
                network["name"],
                RUN_NUM,
                t_num + 1,
                PATTERN_LENGTH,
            )
        )
        with open(filename, "w") as f:
            writer = csv.writer(f, delimiter=",")
            for row in all_pos_sequences[t_num]:
                writer.writerow(row[:pattern_limit])

        filename = (
            DATA_PATH
            + "patterns-activation-peak-%s-run-%03d-trial-%03d-M-%d.csv"
            % (
                network["name"],
                RUN_NUM,
                t_num + 1,
                PATTERN_LENGTH,
            )
        )
        with open(filename, "w") as f:
            writer = csv.writer(f, delimiter=",")
            for row in activation_peak[t_num]:
                writer.writerow(row[:pattern_limit])

        filename = (
            DATA_PATH
            + "patterns-activation-auc-%s-run-%03d-trial-%03d-M-%d.csv"
            % (
                network["name"],
                RUN_NUM,
                t_num + 1,
                PATTERN_LENGTH,
            )
        )
        with open(filename, "w") as f:
            writer = csv.writer(f, delimiter=",")
            for row in activation_auc[t_num]:
                writer.writerow(row[:pattern_limit])

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
