# Wave Response Simulation

The code in this repository simulates the wave response of a multi-modular floating structure and optimizes the stiffness in the connectors. It calculates the system's surge response, relative motion between modules, and the forces in the connecting springs under different wave conditions and spring stiffness values. The results are visualized through various plots.

## Table of Contents
1. [Requirements](#requirements)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Code Structure](#code-structure)
    - [Main Functions](#main-functions)
    - [Helper Functions](#helper-functions)
    - [Additional Files](#additional-files)
5. [Stiffness Optimization](#stiffness-optimization)


## Requirements

The following Python packages are required to run the simulations:

- numpy
- matplotlib
- scipy
- tqdm

## Installation

1. Clone this repository to your local machine.
2. Ensure you have the required Python packages installed (see Requirements).
3. Place the `WadamReader.py` and `waveForce.py` modules in the same directory as the main script.

## Usage

To run the simulations and generate plots, you can use the provided functions in the script.

### Running a Single Simulation

To run a simulation with specific wave period, wave steepness, and spring stiffness values, use the `plotTimeSeriers` function:

## Code Structure

### Main Functions

- `main(WavePeriod, WaveSteepness, springVector, n)`: Runs the main simulation for the specified wave period, wave steepness, and spring stiffness values, where `n` is the number of modules.
- `plotTimeSeriers(WavePeriod, WaveSteepness, springVector)`: Plots the time series response for the specified wave period, wave steepness, and spring stiffness values.
- `runMultipleStiffnesses(modules, stiffnessVector, WavePeriod, WaveSteepness)`: Runs multiple simulations over a range of spring stiffness values and plots the results.

### Helper Functions

- `systemMatrices(n, m_inclAm, d)`: Generates the system matrices `A0`, `B`, `C`, and `Mk` for the simulation.
- `findEigenVal(Mk, stiffnvect, m_inclAm)`: Computes the eigenvalues and eigenfrequencies of the system.

### Additional Files

- `WadamReader.py`: Contains the function `readWadam` for reading WADAM output files.
- `waveForce.py`: Contains the function `w` for calculating the wave force vector in a time step.


# Stiffness Optimization

## Overview

This section provides information about the MATLAB script `alternating2.m` used for optimizing the stiffness in the system based on the simulation code above. The optimization algorithm uses the CVX toolbox for convex optimization in MATLAB.

## Prerequisites
- CVX toolbox installed in MATLAB. You can download and install CVX from [http://cvxr.com/cvx/download/](http://cvxr.com/cvx/download/).

## Installation

1. Ensure the `optimization.m` script and the `SystMatr.mat` file are in the same directory.

## Usage

### Optimization Process

The script performs the following steps:
1. Loads the system matrices from `SystMatr.mat`.
2. Initializes parameters and sets up the optimization problem.
3. Uses CVX to solve a series of semidefinite programming (SDP) problems to find the optimal parameters.
4. Refines the solution iteratively to achieve a tighter bound on the system's performance metric.

## CVX Toolbox

CVX is a MATLAB-based modeling system for convex optimization problems. It allows for easy specification and solution of convex programs. In this script, CVX is used to solve SDP problems during the optimization process.

## Files

- `alternating2.m`: MATLAB script to perform the optimization.
- `SystMatr.mat`: File containing the system matrices `Mk`, `A0`, `B`, and `C`.
