# README

This repository contains minimal code for running the adaptive hybrid rational-LASG algorithm from the paper
>P. Huwiler, D. Pradovera, and J. Schiffmann, _Plug-and-play adaptive surrogate modeling of parametric nonlinear dynamical systems in frequency domain_, International Journal for Numerical Methods in Engineering, 2023

Published version [here](https://doi.org/10.1002/nme.7487). Preprint publicly available [here](http://infoscience.epfl.ch/record/307613?ln=en).

## Disclaimer
The results in the above paper have been obtained by applying the code in this repository to a frequency-domain gas bearing model, which is *not* included in this repository. For more details on the bearing model, we refer to
>E. Guenat and J. Schiffmann, _Effects of humid air on aerodynamic journal bearings_, Tribology International, 2018

Published version [here](https://doi.org/10.1016/j.triboint.2018.06.002). Preprint publicly available [here](http://infoscience.epfl.ch/record/256240?ln=en).

## Prerequisites
* **MATLAB**&reg; (coded and tested on version R2023a)

## Execution
Sample scripts for surrogate training and testing are in `tutorial_train_*.m` and `tutorial_test_*.m`, respectively. The suffixes `1p1` and `1p2` refer to the number of parameters, 1 and 2 (in addition to frequency), respectively.

## Customization
The target QoIs are defined in the function scripts `getSampleImpedance_*.m` (QoIs depending on both frequency and parameters) and `getSampleForce_*.m` (QoIs depending only on parameters). They can be customized to compute any desired QoI. The frequency and parameter ranges are defined near the beginning of `tutorial_train_*.m`, and can be easily adapted to new settings.
