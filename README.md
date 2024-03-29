Reproducing the results of ‘An approximate Bayesian approach for
estimation of the instantaneous reproduction number under misrepoted
epidemic data’
================
Oswaldo Gressani
2022-12-08

[![DOI](https://zenodo.org/badge/532182908.svg)](https://zenodo.org/badge/latestdoi/532182908)

### Reproducibility instructions

This GitHub repository contains all the routines required to reproduce
the tables and figures of the paper entitled “An approximate Bayesian
approach for estimation of the instantaneous reproduction number under
misreported epidemic data”. The repository is structured into four main
folders (**01-Simulations**, **02-Figures**, **03-RealData** and
**04-Complement**). Table 1 below indicates how to reproduce an object
of interest (Object), the location of that object in the main manuscript
(Section), the folder in which to find the desired routine(s) (Folder)
and the R routines (Code(s)) required to reproduce the object of
interest. <br> <br> The routine *simepi.R* is used as a subroutine to
simulate epidemic data and *simul.R* is the core function to reproduce
the simulation results in Section 3.1 (Scenario 1 – Scenario 6) of the
manuscript. The different simulation scenarios have an associated seed,
so that the results can be recovered exactly. In addition, the
**02-Figures** folder contains routines to reproduce Figure 1, 4 and 5,
respectively. Figure 2 and Figure 3 can be recovered with the routines
in the **01-Simulations** folder. The **03-RealData** folder contains
the algorithms to reproduce the real data applications in Section 4 of
the paper (together with Figure 6 and Figure 7). Table 1 has no
associated folder as it contains ‘handwritten’ text that summarizes the
setting of each scenario. Note that the computation time (real elapsed
time) reported in the main manuscript is specific to an Intel Xeon
E-2186M CPU running at a base clock speed of 2.90 GHz.

![](FolderStructure.png)

### Repository version

This is version 0.1.0 (2022-12-08) - “Final structure”.

### License

Copyright © Oswaldo Gressani. All rights reserved.

### URL commit hash

<https://github.com/oswaldogressani/RLPS_misreported_data/commit/4b58443175834a8c40f361c2a86dc74acea8be47>

### Acknowledgments

This project is funded by the European Union’s Research and Innovation
Action under the H2020 work programme, EpiPose (grant number 101003688).
