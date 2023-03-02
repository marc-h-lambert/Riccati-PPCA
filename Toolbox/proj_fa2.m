############################################################################################
###   Companion code for the article 
###   "Low-rank plus diagonal approximations for Riccati-like matrix differential equations"
###   Authors: Silvere Bonnabel, Marc Lambert and Francis Bach
###   Code supported by Silvere Bonnabel and Marc Lambert  
#############################################################################################

# Vectorized variant : shorter but non-optimized for memory
# (Use proj_fa for a memory free algorithm) 
# Project a matrix H=ZZ' Z:d x r on URU'+Diag(barpsi)
function [dU,dR,dbarpsi]  = proj_fa2(Z,U,R,barpsi);
  [d,p]=size(U);
  dbarpsi=proj_diagFA(Z,U);
  
  # Non-optimized computation of diagonal
  M=eye(d)-U*U'
  K=pinv(M.*M);
  mhm=diag(M*(Z*Z')*M);
  dbarpsi=K*mhm;

  G=Z'*U;
  dR=G'*G; # low-rank part
  J=fast_diagU(barpsi,U);
  dR=dR-U'*J; # fa part
  iR=inv(R);
  dU=Z*G*iR-U*G'*G*iR ; # low-rank part
  UPsiU=U'*J;
  dU=dU-J*iR+U*UPsiU*iR; # fa part
end
