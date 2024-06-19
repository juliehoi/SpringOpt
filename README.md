# Wave Response Simulation

The code in this repository simulates the wave response of a multi-modular floating structure. It calculates the system's surge response, relative motion between modules, and the forces in the connecting springs under different wave conditions and spring stiffness values. The results are visualized through various plots.

## Table of Contents
1. [Requirements](#requirements)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Code Structure](#code-structure)
    - [Main Functions](#main-functions)
    - [Helper Functions](#helper-functions)
    - [Additional Files](#additional-files)
5. [Examples](#examples)


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
