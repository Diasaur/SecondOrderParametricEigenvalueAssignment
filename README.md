# SOPEA - SecondOrderParametricEigenvalueAssignment
This repository provides the code for parametric eigenvalue assignment for second order systems.
(Zenodo to this repository: https://doi.org/10.5281/zenodo.15211961)
Originally presented in
https://doi.org/10.1016/j.ymssp.2025.113372
Newest version is submitted for revision

The "matlab_code" folder contains matlab code. Therein
* "Example.mlx" is a live script that is used to run our examples
* "SOPEA.m" contains the algorithm as a function that calls the subcases.
  * "SOPEA_KimuraMinus1.m" contains the new version of the algorithm for m+p=2n.
  * "SOPEA_Kimura.m" contains the new version of the algorithm for m+p=2n+1.
* "SOPEA_Gain.m" calculates the frobenius norm of the resulting gain and is used for optimization.
* "SOPEA_Imag.m" calculates the norm of the imaginary part of the resulting gain and is used for optimization.
* "SecOrd2LTI.m" calculates an LTI system to use matlabs inbuilt eig function
* "SecOrd2Descriptor.m" calculates a Descriptor system to use matlabs inbuilt eig function in case of singular A_2.
* "Feedback2Feedthrough.m" was used to calculate feedback in case of non-zero feedthrough and is now redundant.
