############################################################################################
###   Companion code for the article 
###   "Low-rank plus diagonal approximations for Riccati-like matrix differential equations"
###   Authors: Silvere Bonnabel, Marc Lambert and Francis Bach
###   Code supported by Silvere Bonnabel and Marc Lambert  
#############################################################################################

#############################################################################################
# Brute force Riccati with FA decomposition P=URU' + Psi
# Attention: Simpler to read but not memory efficient !! Use Riccati_fa_fast instead
#############################################################################################
# U low rank matrix d x p
# R matrix p x p (diagonal)
# diagPsi vector of dim d correponding to the diagonal of Psi
# dt integration step 
# sqrtA sqrt of A st sqrtA*sqrtA'=A and sqrtA d x k (k<<d)
# --> if A not symmetric use a sparse structure 
# C obs matrix m x d
# diagQ diagonal of the process noise matrix d x d 
# N obs noise matrix m x m
#############################################################################################
 
function [U,R,diagPsi] = Riccati_fa(U,R,diagPsi,dt,sqrtA,C,dq,N);
  Psi=diag(diagPsi);
  Q=diag(dq);
  A=sqrtA*sqrtA';
  [d,p]=size(U);
  I=eye(d);
  M=I-U*U';
  S=C'*inv(N)*C;
  M1=A*Psi+Psi*A'+Q-Psi*S*Psi;
  M2=R*U'*S*Psi*U;
  M3=U'*A*U*R;
    
  #compute dPsi
  K=pinv(M.*M);
  mhm=diag(M*(M1)*M);
  dbarPsi=K*mhm;
  dPsi=diag(dbarPsi);
    
  dU=M*(M1-dPsi)*U*inv(R)+M*(A*U-Psi*S*U); 
  dR=U'*(M1-dPsi)*U-M2-M2'-R*U'*S*U*R+M3+M3';
    
  # project on manifold
  U=retraction_qr(U,dU,dt);
  R=R+dt*(dR);
  R=(R+R')/2;
  Psi=Psi+dt*dPsi;
  diagPsi=diag(Psi);
end