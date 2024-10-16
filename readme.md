This readme describes the steps to reproduce the paper figures


## Software versions
- PYHTON : 3.10.12 
- NEURON : 8.2.3 
- MOOSE : 4.0.0 (pymoose)
- Version numbers of all python packages are listed in the requirements.txt file present in the repository. 

## Setting up
Git clone the repository and go into synaptic-convergence-paper directory and execute the below lines
```
python3 -m venv conv-env
source conv-env/bin/activate
pip install -r requirements.txt
```


## Figure 1

Only schematics

### Plot generation

Path : synaptic-convergence-paper/simulations/figures/figure1\
command : ```python3 plot_figure1.py```\
_Ensure that schematics are present in the directory containing the plotting file_



## Figure 2

### Data generation

- Figure 2A, C, D, E\
  Schematics\
  _Schematics need to be present in the directory containing the plotting file_
- Figure 2F, G, H, I, J, K, M\
  Path : synaptic-convergence-paper/simulations/codes/groups\
  command : ```python3 group_counts.py <run_number> <seed>```\
  _The run number and seeds can be found in the same path. It takes around ~20 hrs for each run-number and seed combination. 100 processes were launched in parallel with different seeds for the analysis mentioned in the paper._
- Figure 2B\
  Path :  synaptic-convergence-paper/simulations/codes/groups/cam\
  command : ```python3 generate_Ca_CaM_trace.py```\
  _Ensure that the model file Ca_CaM.g is present in the same path_
- Figure 2L\
  Path : synaptic-convergence-paper/simulations/codes/groups\
  command : ```python3 categorize_stimulus_driven_groups.py <run_number> <seed>```\
  _The run number and seeds can be found in the same path. It takes around 4 mins for each run-number and seed combination. 100 processes were launched in parallel with different seeds for the analysis mentioned in the paper._

### Plot generation

Path : synaptic-convergence-paper/simulations/figures/figure2\
command : ```python3 plot_figure2.py```


## Figure 2 - Supplement 1

### Plot generation

Path : synaptic-convergence-paper/simulations/figures/figure2\
command : ```python3 plot_figure2_param_screen.py```



## Figure 3

### Data generation

- Figure 3B, C, D, E, F, G, H, I, J, K\
  Path : synaptic-convergence-paper/simulations/codes/groups\
  command : ```python3 group_counts_with_activation.py <run_number> <seed>```\
  _The run number and seeds can be found in the same path. It takes around 5hrs for each run-number and seed combination for hippo-elec. 100 processes were launched in parallel with different seeds for the analysis mentioned in the paper. Might take longer for some of the other network configurations_

### Plot generation

Path : synaptic-convergence-paper/simulations/figures/figure3\
command : ```python3 plot_figure3.py```



## Figure 4

### Data generation

- Figure 4A, B, C\
  Schematics\
  _Schematics needs to be present in the directory containing the plotting file_
- Figure 4D, E, F, G, H, I\
  Path : synaptic-convergence-paper/simulations/codes/sequences\
  command: ```python3 sequence_counts.py <run_number> <seed>```\
  _The run number and seeds can be found in the same path. It takes around ~4 hrs for each run-number and seed combination. 100 processes were launched in parallel for the analysis mentioned in the paper._

### Plot generation

Path : synaptic-convergence-paper/simulations/figures/figure4\
command : ```python3 plot_figure4.py```


## Figure 4 - Supplement 1

### Plot generation

Path : synaptic-convergence-paper/simulations/figures/figure4\
command : ```python3 plot_figure4_param_screen.py```



## Figure 5

### Data generation

- Figure 5A\
  Schematic\
  _Schematic needs to be present in the directory containing the plotting file_
- Figure 5B\
  Path : synaptic-convergence-paper/simulations/codes/sequences/chemical\
  command : ```python3 sim_seq_vs_scrambled.py```
- Figure 5C\
  Path : synaptic-convergence-paper/simulations/codes/sequences/chemical\
  command : ```python3 sim_q_scores_vs_stimulus_patterns.py```
- Figure 5D\
  Schematic\
  _Schematic needs to be present in the directory containing the plotting file_
- Figure 5E\
  Path : synaptic-convergence-paper/simulations/codes/sequences/chemical\
  command : ```python3 sim_ectopic_positions.py```\
  runtime : 70 min
- Figure 5F\
  Path : synaptic-convergence-paper/simulations/codes/sequences/chemical\
  command : ```python3 sim_selectivity_seq_length.py```

### Plot generation

Path : synaptic-convergence-paper/simulations/figures/figure5\
command : ```python3 plot_figure5.py```



## Figure 6

__Note__\
_For the below figure we chose the model published on ModelDB (https://modeldb.science/140828) for examining the effect of ectopic inputs on the selectivity of electrical sequences (Branco, Clark, and Häusser 2010). The base code from ModelDB was modified and run using the NEURON Simulator. The original channel kinetics parameters post the fix to NMDA_Mg_T initialization, and the dendrite used for testing were retained as is from ModelDB github repository (https://github.com/ModelDBRepository/140828). We tested out the effects of ectopics using the passive version of the model._

### Data generation

- Figure 6A\
  Schematic\
  _Schematic needs to be present in the directory containing the plotting file_
- Figure 6B\
  Path : synaptic-convergence-paper/simulations/codes/sequences/electrical\
  command : ```python3 main.py```\
  _Ensure that all the dependent files are present. Sequence of length 5 was chosen for this._
- Figure 6C\
  Path : synaptic-convergence-paper/simulations/codes/sequences/electrical\
  command : ```python3 main.py```\
  _The command needs to be run separately for sequences of 3, 5 and 8 inputs. Ensure that all the dependent files for each sequence length are present. Accordingly modify the synapse_loc_<sequence_length>.dat filename in the init_synapses.hoc file and the corresponding synapse_order_<sequence_length>.dat filename in init_sim.hoc. The output needs to be redirected to tempV_<sequence_length>_panelC.dat and tempM_<sequence_length>_panelC.dat files as well in init_sim.hoc._
- Figure 6D\
  Path : synaptic-convergence-paper/simulations/codes/sequences/electrical\
  command : ```python3 main.py```\
  runtime : ~50 min
  _Note that the mosinit.hoc file will need to call init_sim_ectopic.hoc instead of init_sim.hoc. This was done for a sequence of length 5. Accordingly use synapse_order_5.dat, ectopic_input_time_panelD.dat in init_sim_ectopic.hoc. Use synapse_loc_5.dat and ectopic_synapse_loc_panelD.dat in init_synapses_ectopic.hoc. Use The output needs to be redirected to tempV_panelD.dat and tempM_panelD.dat as well in init_sim_ectopic.hoc._
- Figure 6E\
  Path : synaptic-convergence-paper/simulations/codes/sequences/electrical\
  command : ```python3 main.py```\
  runtime : ~30 min in total for sequences of length 3 to 9
  _Note that the mosinit.hoc file will need to call init_sim_ectopic.hoc instead of init_sim.hoc. The command needs to be run separately for sequences of different sequence lengths. Ensure that all the dependent files for each sequence length are present. Accordingly use synapse_order_<sequence_length>.dat, ectopic_input_time_<sequence_length>.dat in init_sim_ectopic.hoc. Use synapse_loc_<sequence_length>.dat and ectopic_synapse_loc_<sequence_length>.dat in init_synapses_ectopic.hoc. Use The output needs to be redirected to tempV_<sequence_length>_panelE.dat and tempM_<sequence_length>_panelE.dat as well in init_sim_ectopic.hoc._

### Plot generation

Path : synaptic-convergence-paper/simulations/figures/figure6\
command : ```python3 plot_figure6.py```


## Figure 6 - Supplement 1

### Data generation

_Follow instructions in Figure 6D_

### Plot generation

Path : synaptic-convergence-paper/simulations/figures/figure6\
command : ```python3 plot_ectopic_responses.py```



## Figure 7

### Data generation

- Figure 7 - all panels\
  Path : synaptic-convergence-paper/simulations/codes/sequences\
  command : ```python3 seq_selectivity_with_activation.py <run_number> <seed>```\
  _The run number and seeds can be found in the same path._\
  _This takes around ~13 hrs when for each run-number and seed combination for hippo-CICR for T_POP_POST=10000 neurons. 100 processes were launched in parallel for the analysis mentioned in the paper._

### Plot generation

Path : synaptic-convergence-paper/simulations/figures/figure7\
command : ```python3 plot_figure7.py```


## Figure 7 - Supplement 1

### Data generation

_Follow instructions in Figure 7. In seq_selectivity_with_activation.py, use prm.cortex_cicr in networks and T_POP_POST=4000_

### Plot generation

Path : synaptic-convergence-paper/simulations/figures/figure7\
command : ```python3 plot_figure7_cortex_cicr.py```


## Figure 8

**Author** - Upinder S Bhalla

### Data generation

- Figure 8B, C, D, E, F, G\
  Path : synaptic-convergence-paper/simulations/figures/figure8\
  command: 
  ```
  ./runMAPKSeries.bat
  ./runVmSeries.bat  
  ./runCaMSeries.bat  
  ```

### Plot generation

Path : synaptic-convergence-paper/simulations/figures/figure8\
command : ```python3 analyzeNonlin6.py```



