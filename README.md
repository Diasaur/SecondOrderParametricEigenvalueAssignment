# SOPEA - SecondOrderParametricEigenvalueAssignment
This repository provides the code for parametric eigenvalue assignment for second order systems.

The "matlab_code" folder contains matlab code.
* Therein the "Example.mlx" is a live script that is used to run our examples
* "ParametricEigenvalueAssignment.m" contains the algorithm as a function
* "Feedback2Feedthrough.m" is used to calculate feedback in case of non-zero feedthrough
* "ParametricEigenvalueAssignmentGain.m" calculates the frobenius gain and is used for optimization
* "SecOrd2LTI.m" calculates an LTI system to use matlabs inbuilt eig function
