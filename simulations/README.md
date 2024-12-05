# mixWAS Simulations
Simulations for mixWAS. All simulations were run on the FAS Research Computing cluster at Harvard University using R 4.2.2.

## Simulation Scripts
* __generate_data.R__: Function to generate data for mixWAS simulations
* __run_simulation.R__: Wrapper to conduct the simulation on a given list of input (e.g. fixed $\beta$)
* __pheWAS.R__: Contains functions for running PheWAS mega/meta
* __helpers.R__: Helper functions
* __power_simulation.R__: Wrapper to run simulations for all parameters (e.g. range of $\beta$), on the cluster
* __simulation_inputs_v3.R__: Function to define all simulation inputs.

## Inputs/Results
* __inputs/__: Folder of simulation inputs
* __outputs/__: Folder of results
* __paper_figues.R__: Script to generate simulation figures for paper.
* __paper_figures/__: Simulation figures in paper

## Mapping of Simulation Inputs to Paper Simulation

Mapping of simulation IDs to their corresponding ordering in the final version of the paper.

| Simulation ID |                            Paper Simulation Description                           |
|:-------------:|:---------------------------------------------------------------------------------:|
|  same_sign_1  |      Simulation 1: $\beta$ in same direction w/ positive residual correlation     |
|  same_sign_2  |         Simulation 1: $\beta$ in same direction w/ no residual correlation        |
|  same_sign_3  |      Simulation 1: $\beta$ in same direction w/ negative residual correlation     |
|   opp_sign_1  |   Simulation 2: $\beta$ in opposite directions w/ positive residual correlation (SUPPLEMENT)   |
|   opp_sign_2  |      Simulation 2: $\beta$ in opposite directions w/ no residual correlation  (SUPPLEMENT)     |
|   opp_sign_3  |   Simulation 2: $\beta$ in opposite directions w/ negative residual correlation (SUPPLEMENT)   |
|       2       |         Simulation 3: Binary Phenotypes Only, Common Variant (SUPPLEMENT)         |
|       3       |          Simulation 4: Binary Phenotypes Only, Rare Variant (SUPPLEMENT)          |
|       4       | Simulation 5: Binary Phenotypes Only, Healthy Controls (SUPPLEMENT) |

