# Riccati-PPCA
Companion code for the article "Low-rank plus diagonal approximations for Riccati-like matrix differential equations". Authors: Silvere Bonnabel, Marc Lambert and Francis Bach.

These programs compute the low-rank $URU^T$ and low rank + diagonal approximation $URU^T+\psi$ of a matrix $H$ supposed of the form $H=GG^T$. The PPCA approximation is when $\psi$ is of the form $s(\mathbb{I}-UU^T)$ with $s$ a scalar. The FA approximation is when the matrix $\psi$ is diagonal.

- [XP_Riccati_main][1]: Run to reproduce the figures of the paper. Successive projections to track a time dependent covariance matrix $H=P_t$ given by the Ricatti equations. The sequence of projections are given by ODE on each of the factorization components : $U_t$, $R_t$ and optionally $\psi_t$. 
May takes some time, the observation matrix is of large dimension here.

- [XP_Riccati_fast][1]: Experiment on a simpler exemple with an observation matrix with a limited number of rows.

- [XP_Proj][3]: Experiment on a single step projection in dimension one million.

- [FactorAnalysis][4]: This directory contains the mains methods for FA-PPCA projection and FA-PPCA Riccati ODE. Riccati function are given in  two versions: a non memory-efficient variant easier to read for validation purpose and a fast memory-efficient variant. The low-rank and full methods are also given. 

- [Toolbox][5]: This directory contains usefull function for fast linear lagebra on diagonal and square root matrix and also the retraction function.

[1]: ./XP_Riccati_main.m
[2]: ./AdditionalXP/XP_Riccati_fast.m
[3]: ./AdditionalXP/XP_Proj.m
[4]: ./FactorAnalysis 
[5]: ./Toolbox
