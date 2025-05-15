# Subgradient Regularization Methods for Nonsmooth Optimization

MATLAB code for the numerical experiments in the paper:
**"Subgradient Regularization: A Descent-Oriented Subgradient Method for Nonsmooth Optimization"**
Hanyang Li, Ying Cui
[https://arxiv.org/abs/2505.07143](https://arxiv.org/abs/2505.07143)

## Overview

This repository provides MATLAB implementations for testing nonsmooth optimization algorithms, focusing on our proposed Subgradient-Regularized Descent methods (SRD, SRD-adapt). The code facilitates replication of the numerical results presented in our paper.

Key components:
* **Experiment Scripts:** Located in `experiments/`, these scripts define specific test scenarios (problems, parameters) and initiate runs.
* **Core Engine:** `src/run_single_experiment.m` manages the execution of these scenarios.
* **Solvers:** Algorithm implementations are in `src/core_solvers/`.
* **Problems:** Test function definitions are in `src/problem_definitions/`.
* **Outputs:** Results (data, tables, logs) are saved in the `results/` directory, organized by experiment group.

## Requirements

* MATLAB R2020b or later (tested with R2023a).
* Optimization Toolbox (for `quadprog`, used by some solvers).

## Algorithm Notes & External Code

For comparison, this codebase also includes or adapts the following algorithms. Please refer to their original sources for complete details:

* **Polyak:** Standard subgradient method with Polyak stepsizes.
* **PBMDC:** A proximal bundle method for DC functions. Based on code from [www.oliveira.mat.br/solvers](https://www.oliveira.mat.br/solvers).
* **NTD:** Our MATLAB re-implementation of an NTD method for fair runtime comparison, originally in PyTorch. *(Original source: [github.com/COR-OPT/ntd.py](https://github.com/COR-OPT/ntd.py))*
* **HANSO:** A hybrid nonsmooth optimization algorithm (BFGS and Gradient Sampling). Based on the original from [cs.nyu.edu/~overton/software/hanso/](https://cs.nyu.edu/~overton/software/hanso/). We test its BFGS and GS components separately.
