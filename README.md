# SOPEA - SecondOrderParametricEigenvalueAssignment
This repository provides the code for parametric eigenvalue assignment for second order systems.
https://doi.org/10.5281/zenodo.15211961
As presented in
https://doi.org/10.1016/j.ymssp.2025.113372

The "matlab_code" folder contains matlab code.
* Therein the "Example.mlx" is a live script that is used to run our examples
* "ParametricEigenvalueAssignment.m" contains the algorithm as a function
* "Feedback2Feedthrough.m" is used to calculate feedback in case of non-zero feedthrough
* "ParametricEigenvalueAssignmentGain.m" calculates the frobenius gain and is used for optimization
* "SecOrd2LTI.m" calculates an LTI system to use matlabs inbuilt eig function
