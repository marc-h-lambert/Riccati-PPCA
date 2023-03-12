# Riccati-PPCA
Companion code for the article "Low-rank plus diagonal approximations for Riccati-like matrix differential equations". Authors: Silvere Bonnabel, Marc Lambert and Francis Bach.

These programs compute the low-rank $URU^T$ and low rank + diagonal approximation $URU^T+\psi$ of a matrix $H$ supposed of the form $H=GG^T$. The PPCA approximation is when $\psi$ is of the form $s\mathbb{I}$ with $s$ a scalar. The FA approximation is when the matrix $\psi$ is diagonal.

- [XP_Riccati_comparison][1]: Run to reproduce the figures of the paper. Successive projections to track a time dependent covariance matrix $H=P_t$ given by the Ricatti equations. The sequence of projections are given by ODE on each of the factorization components : $U_t$, $R_t$ and optionally $\psi_t$. 

- [XP_Proj][2]: Single projection of a matrix $H$ on the tangent plane of a low-rank manifold. (Used to generate the Table 1 in the paper)

[1]: ./XP_Riccati_comparison.m
[2]: ./XP_Proj.m
[3]: ./XP_Riccati.m
