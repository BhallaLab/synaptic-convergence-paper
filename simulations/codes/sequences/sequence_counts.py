import numpy as np
import csv
import time
import sys
import sequence_params as prm
import os
from numba import njit
from scipy import sparse

np.set_printoptions(threshold=np.inf)


MIN_M = prm.MIN_M
MAX_M = prm.MAX_M

T_POP_POST = 4000

RUN_NUM = int(sys.argv[1])
SEED = int(sys.argv[2])
STIMULUS_SEED = 47
NUM_TRIALS = 100
STIM_ORDER = np.arange(1, MAX_M + 1)

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
    stim_order,
    ensemble_size,
    ensemble_participation_prob,
    stim_type,
):
    """Simulate presynaptic activity for the duration of a MAX_M sequence time steps

    Parameters
    ----------
    t_pop_pre : int
        Size of the presynaptic population
    noise_prob : float
        Probability of background activity
    stim_order : array like
        Order of stimulation of the ensembles
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
        2D array of the size of (t_pop_pre x MAX_M), representing presynaptic activity.
        -1: Active due to background activity
        0: Inactive
        1: Active due to participation in the ensemble
    """
    match stim_type:
        case 0:
            pre_pop_activity = np.random.choice(
                np.arange(-1, 1, dtype=np.int8),
                size=(t_pop_pre, MAX_M),
                p=[noise_prob, 1 - noise_prob],
            )
        case 1:
            pre_pop_activity = np.random.choice(
                np.arange(-1, 1, dtype=np.int8),
                size=(t_pop_pre, MAX_M),
                p=[noise_prob, 1 - noise_prob],
            )
            for i, ensemble_id in enumerate(stim_order):
                pre_pop_activity[
                    (ensemble_id - 1) * ensemble_size : ensemble_id * ensemble_size, i
                ] = np.random.choice(
                    [1, 0],
                    size=ensemble_size,
                    p=[ensemble_participation_prob, 1 - ensemble_participation_prob],
                )
        case -1:
            pre_pop_activity = np.zeros((t_pop_pre, MAX_M), dtype=np.int8)
            for i, ensemble_id in enumerate(stim_order):
                pre_pop_activity[
                    (ensemble_id - 1) * ensemble_size : ensemble_id * ensemble_size, i
                ] = 1
        case _:
            pre_pop_activity = np.zeros((t_pop_pre, MAX_M), dtype=np.int8)
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


def count_sequences(t_pop_pre, N, p_conn, p_e, L, S, delta, R, D):
    """Count different kinds of sequences occurring on a neuron's dendrite

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

    Returns
    -------
    stimulus : np.ndarray
        1D-array of stimulus type of size number of trials containing values from [-1, 0, 1].
        -1: All ensemble neurons are active due to the stimulus
        0: No ensemble activity, only background activity
        1: Ensembles with a fraction of neurons active, along with background activity in the non-ensemblar neurons
    poss_sequences : np.ndarray
        Array of size (num_trials x postsyn_population_size x num_seq_lengths) containing the number of poss sequences on a neuron
    noise_sequences : np.ndarray
        Array of size (num_trials x postsyn_population_size x num_seq_lengths) containing the number of noise sequences on a neuron
    all_sequences : np.ndarray
        Array of size (num_trials x postsyn_population_size x num_seq_lengths) containing the number of all sequences on a neuron
    """
    num_synapses = int(p_conn * t_pop_pre)
    spacing_size = int(S * num_synapses / L)
    input_width = int(delta * num_synapses / L)
    noise_prob = 1 - np.exp(-R * D)
    poss_sequences = np.zeros((NUM_TRIALS, T_POP_POST, MAX_M + 1), dtype=np.uint64)
    noise_sequences = np.zeros((NUM_TRIALS, T_POP_POST, MAX_M + 1), dtype=np.uint64)
    all_sequences = np.zeros((NUM_TRIALS, T_POP_POST, MAX_M + 1), dtype=np.uint64)

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
            t_pop_pre, noise_prob, STIM_ORDER, N, p_e, stimulus[t_num]
        )
        for neuron_id in range(T_POP_POST):
            curr_dist_matrix_noise_seq = sparse.csr_array(
                gen_dist_matrix(
                    "noise",
                    activity[:, 0],
                    activity[:, 1],
                    pre_conn[neuron_id],
                    num_synapses,
                    spacing_size,
                    input_width,
                ),
                shape=(num_synapses, num_synapses),
            )
            cum_paths_noise_seq = curr_dist_matrix_noise_seq
            noise_sequences[t_num, neuron_id, 2] = np.sum(cum_paths_noise_seq)

            if stimulus[t_num]:
                curr_dist_matrix_poss_seq = sparse.csr_array(
                    gen_dist_matrix(
                        "poss",
                        activity[:, 0],
                        activity[:, 1],
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
                        activity[:, 0],
                        activity[:, 1],
                        pre_conn[neuron_id],
                        num_synapses,
                        spacing_size,
                        input_width,
                    ),
                    shape=(num_synapses, num_synapses),
                )
                cum_paths_all_seq = curr_dist_matrix_all_seq
                poss_sequences[t_num, neuron_id, 2] = np.sum(cum_paths_poss_seq)
                all_sequences[t_num, neuron_id, 2] = np.sum(cum_paths_all_seq)

            for M in range(3, MAX_M + 1):
                if noise_sequences[t_num, neuron_id, M - 1] > 0:
                    next_dist_matrix_noise_seq = sparse.csr_array(
                        gen_dist_matrix(
                            "noise",
                            activity[:, M - 2],
                            activity[:, M - 1],
                            pre_conn[neuron_id],
                            num_synapses,
                            spacing_size,
                            input_width,
                        ),
                        shape=(num_synapses, num_synapses),
                    )
                    cum_paths_noise_seq = mul_mat(
                        cum_paths_noise_seq, next_dist_matrix_noise_seq
                    )
                    noise_sequences[t_num, neuron_id, M] = np.sum(cum_paths_noise_seq)
                    curr_dist_matrix_noise_seq = next_dist_matrix_noise_seq

                if stimulus[t_num]:
                    if poss_sequences[t_num, neuron_id, M - 1] > 0:
                        next_dist_matrix_poss_seq = sparse.csr_array(
                            gen_dist_matrix(
                                "poss",
                                activity[:, M - 2],
                                activity[:, M - 1],
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
                        poss_sequences[t_num, neuron_id, M] = np.sum(cum_paths_poss_seq)
                        curr_dist_matrix_poss_seq = next_dist_matrix_poss_seq

                    if all_sequences[t_num, neuron_id, M - 1] > 0:
                        next_dist_matrix_all_seq = sparse.csr_array(
                            gen_dist_matrix(
                                "all",
                                activity[:, M - 2],
                                activity[:, M - 1],
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
                        all_sequences[t_num, neuron_id, M] = np.sum(cum_paths_all_seq)
                        curr_dist_matrix_all_seq = next_dist_matrix_all_seq

    return stimulus, poss_sequences, noise_sequences, all_sequences


def process_network(network):
    """Simulate presynaptic activity and count sequences for the network

    Parameters
    ----------
    network : dict
        Dictionary containing the parameters of the network
    """
    (
        stimulus,
        poss_pos_sequences,
        noise_pos_sequences,
        all_pos_sequences,
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
    )
    filename = DATA_PATH + "stimulus-%s-run-%03d.csv" % (network["name"], RUN_NUM)
    with open(filename, "w") as f:
        writer = csv.writer(f, delimiter=",")
        writer.writerow(stimulus)
    for t_num in range(NUM_TRIALS):
        filename = DATA_PATH + "noise-sequences-%s-run-%03d-trial-%03d.csv" % (
            network["name"],
            RUN_NUM,
            t_num + 1,
        )
        with open(filename, "w") as f:
            writer = csv.writer(f, delimiter=",")
            for row in noise_pos_sequences[t_num]:
                writer.writerow(row[MIN_M:])

        if stimulus[t_num]:
            filename = (
                DATA_PATH
                + "perfect-stimulus-driven-sequences-%s-run-%03d-trial-%03d.csv"
                % (
                    network["name"],
                    RUN_NUM,
                    t_num + 1,
                )
            )
            with open(filename, "w") as f:
                writer = csv.writer(f, delimiter=",")
                for row in poss_pos_sequences[t_num]:
                    writer.writerow(row[MIN_M:])

            filename = DATA_PATH + "all-sequences-%s-run-%03d-trial-%03d.csv" % (
                network["name"],
                RUN_NUM,
                t_num + 1,
            )
            with open(filename, "w") as f:
                writer = csv.writer(f, delimiter=",")
                for row in all_pos_sequences[t_num]:
                    writer.writerow(row[MIN_M:])

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
