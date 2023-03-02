# Riccati-PPCA
Companion code for the article "Low-rank plus diagonal approximations for Riccati-like matrix differential equations". Authors: Silvere Bonnabel, Marc Lambert and Francis Bach.

These programs compute the low-rank URU' and low rank + diagonal approximation URU'+Psi of a matrix H supposed of the form H=GG'. The PPCA approximation is when Psi is of the form s.Id with s a scalar. The FA approximation is when the matrix is diagonal.

XP_Proj: Single projection of a matrix H on the tangent plane of a low-rank manifold.

XP_Ricatti: Successive projections to track a time dependent covariance matrix P given by the Ricatti equations. The sequence of projections are given by ODE on each of the factorization components : U, R and optionally Psi. 
