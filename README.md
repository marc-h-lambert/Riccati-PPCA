# Riccati-PPCA
Companion code for the article "Low-rank plus diagonal approximations for Riccati-like matrix differential equations". Authors: Silvere Bonnabel, Marc Lambert and Francis Bach.

These programs compute the projections of a structured symmetric matrix $H$ onto the tangent spaces to the manifolds of low-rank $URU^T$ and low rank + diagonal approximation $URU^T+\psi$ (FA) and low-rank + isotropic diagonal $URU^T+s(\mathbb{I}-UU^T)$ (PPCA) matrices.  

- [XP_Proj][1]: Experiment on a single step projection in dimension one million of a structured symmetric matrix $H=GG^T$.

- [XP_Riccati_main][2]: Run to reproduce the figures of the paper. Successive projections to track a time dependent covariance matrix $H=P_t$ given by the Ricatti equations. The sequence of projections are given by ODE on each of the factorization components : $U_t$, $R_t$ and optionally $\psi_t$.

- [FactorAnalysis][3]: CHANGER LE NOM ? This directory contains the mains methods for FA and PPCA projection and FA and PPCA Riccati ODE. Riccati function are given in  two versions: a non memory-efficient variant easier to read for validation purpose and a fast memory-efficient variant. The low-rank and full methods are also given.

- [Toolbox][4]: This directory contains usefull function for fast linear lagebra on diagonal and square root matrix and also the retraction function.

[1]: ./XP_Riccati_main.m
[2]: ./XP_Riccati_fast.m
[3]: ./XP_Proj.m
[4]: ./FactorAnalysis 
[5]: ./Toolbox
