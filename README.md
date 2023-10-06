# README

This repository contains minimal code for running the adaptive hybrid rational-LASG algorithm from the upcoming paper:

P. Huwiler, D. Pradovera, and J. Schiffmann, _Plug-and-play adaptive surrogate modeling of parametric nonlinear dynamical systems in frequency domain_ (2023)

## Prerequisites
* **MATLAB**&reg; (coded and tested on version R2023a)

## Execution
Sample scripts for surrogate training and testing are in `tutorial_train_*.m` and `tutorial_test_*.m`, respectively. The suffixes `1p1` and `1p2` refer to the number of parameters, 1 and 2 (in addition to frequency), respectively.

## Customization
The target QoIs are defined in the function scripts `getSampleImpedance.m` (QoIs depending on both frequency and parameters) and `getSampleForce.m` (QoIs depending only on parameters). They can be customized to compute any desired QoI. The frequency and parameter ranges are defined near the beginning of `tutorial_train_*.m`, and can be easily adapted to new settings.
