############################################################################################
###   Companion code for the article 
###   "Low-rank plus diagonal approximations for Riccati-like matrix differential equations"
###   Authors: Silvere Bonnabel, Marc Lambert and Francis Bach
###   Code supported by Silvere Bonnabel and Marc Lambert  
#############################################################################################

#############################################################################################
# Brute force Riccati with low-rank decomposition P=URU'
# Attention: Simpler to read but not memory efficient !! Use Riccati_fa_fast instead
#############################################################################################
# U low rank matrix d x p
# R matrix p x p (diagonal)
# dt integration step 
# sqrtA sqrt of A st sqrtA*sqrtA'=A and sqrtA d x k (k<<d)
# --> if A not symmetric use a sparse structure 
# C obs matrix m x d
# diagQ diagonal of the process noise matrix d x d 
# N obs noise matrix m x m
#############################################################################################
 
 
function [U,R] = Riccati_lowRank(U,R,dt,sqrtA,C,dq,N);
  Q=diag(dq);
  A=sqrtA*sqrtA';
  [d,p]=size(U);
  dU=(eye(d)-U*U')*(A*U+Q*U*inv(R));
  M=U'*A*U*R;
  dR=M+M'+U'*Q*U-R*U'*C'*inv(N)*C*U*R;
  U=retraction_qr(U,dU,dt);
  R=R+dt*(dR);
  R=(R+R')/2;
end