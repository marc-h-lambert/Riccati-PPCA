############################################################################################
###   Companion code for the article 
###   "Low-rank plus diagonal approximations for Riccati-like matrix differential equations"
###   Authors: Silvere Bonnabel, Marc Lambert and Francis Bach
###   Code supported by Silvere Bonnabel and Marc Lambert  
#############################################################################################

# Project a matrix H=ZZ' Z:d x r on URU'+Diag(barpsi)
function [dU,dR,dbarpsi]  = proj_fa(Z,U,R,barpsi);
  [d,p]=size(U);
  d=fast_diagFA1(Z,Z,U);
  dbarpsi=solve_diagFA(d,U);
  G=Z'*U;
  dR=G'*G; # low-rank part
  J=fast_diagU(barpsi,U);
  dR=dR-U'*J; # fa part
  iR=inv(R);
  dU=Z*G*iR-U*G'*G*iR ; # low-rank part
  UPsiU=U'*J;
  dU=dU-J*iR+U*UPsiU*iR; # fa part
end
