############################################################################################
###   Companion code for the article 
###   "Low-rank plus diagonal approximations for Riccati-like matrix differential equations"
###   Authors: Silvere Bonnabel, Marc Lambert and Francis Bach
###   Code supported by Silvere Bonnabel and Marc Lambert  
#############################################################################################


#############################################################################################
# Brute force Riccati with PPCA decomposition P=URU' + s(I-UU')
# Attention: Simpler to read but not memory efficient !! Use Riccati_fa_fast instead
#############################################################################################
# U low rank matrix d x p
# R matrix p x p (diagonal)
# s scalar driving the isotropic variance term of P
# dt integration step 
# sqrtA sqrt of A st sqrtA*sqrtA'=A and sqrtA d x k (k<<d)
# --> if A not symmetric use a sparse structure 
# C obs matrix m x d
# diagQ diagonal of the process noise matrix d x d 
# N obs noise matrix m x m
#############################################################################################
 
 
function [U,R,s] = Riccati_ppca(U,R,s,dt,sqrtA,C,dq,N);
  Q=diag(dq);
  A=sqrtA*sqrtA';
  [d,p]=size(U);
  S=C'*inv(N)*C;
  dsigma=trace((eye(d)-U*U')*(2*s*A+Q-s*s*S))/(d-p);
  dU=(eye(d)-U*U')*(A*U*R+Q*U+s*A'*U-s*S*U*R)*inv(R-s*eye(p));
  dR=U'*A*U*R+R*U'*A*U+U'*Q*U-R*U'*S*U*R;
  U=retraction_qr(U,dU,dt);
  R=R+dt*(dR);
  R=(R+R')/2;
  s=s+dt*dsigma;
end