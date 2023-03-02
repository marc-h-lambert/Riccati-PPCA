# Riccati-PPCA
Companion code for the article "Low-rank plus diagonal approximations for Riccati-like matrix differential equations". Authors: Silvere Bonnabel, Marc Lambert and Francis Bach.

These programs compute the low-rank URU' and low rank + diagonal approximation $URU^T+\psi$ of a matrix H supposed of the form $H=GG^T$. The PPCA approximation is when $\psi$ is of the form $s\mathbb{I}$ with $s$ a scalar. The FA approximation is when the matrix $\psi$ is diagonal.

XP_Proj: Single projection of a matrix $H$ on the tangent plane of a low-rank manifold.

XP_Riccati: Successive projections to track a time dependent covariance matrix $P$ given by the Ricatti equations. The sequence of projections are given by ODE on each of the factorization components : $U$, $R$ and optionally $\psi$. 
