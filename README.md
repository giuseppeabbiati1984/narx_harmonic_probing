# Numerical Harmonic Probing of Polynomial NARX Models

This repository contains MATLAB code for estimating generalized frequency response functions (GFRFs) from a polynomial NARX model using the numerical harmonic probing algorithm. The example problem is a nonlinear Duffing oscillator, for which synthetic input-output data are generated, a polynomial NARX model is trained, and the resulting model is probed to estimate the linear, quadratic, and cubic transfer functions.

The implementation follows the numerical harmonic probing framework described in:

> Stamenov et al. (2025), *Numerical estimation of generalized frequency response functions from time series data using NARX*, Mechanical Systems and Signal Processing.  
> https://doi.org/10.1016/j.ymssp.2025.113278

## Overview

The code demonstrates the following workflow:

1. Generate synthetic data from a Duffing oscillator.
2. Segment the data for polynomial NARX model training.
3. Train one polynomial NARX model for each data segment.
4. Extract the first-, second-, and third-order NARX coefficient arrays.
5. Apply numerical harmonic probing to estimate:
   - the linear transfer function, `H1`;
   - the quadratic transfer function, `H2`;
   - the cubic transfer function, `H3`.
6. Compare the estimated transfer functions against analytical Duffing oscillator solutions.
7. Plot real and imaginary parts of the estimated and theoretical transfer functions.

## Repository Structure

├── script_data_gen.m
├── script_harmonic_probing.m
├── compute_probing_H1.m
├── compute_probing_H2.m
├── compute_probing_H3.m
├── compute_freq_check.m
├── find_freq.m
├── README.md
└── data/
    ├── duff_train_data.mat
    └── DuffTF_dw=4.mat
