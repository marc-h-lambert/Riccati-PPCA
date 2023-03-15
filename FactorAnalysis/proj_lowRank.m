############################################################################################
###   Companion code for the article 
###   "Low-rank plus diagonal approximations for Riccati-like matrix differential equations"
###   Authors: Silvere Bonnabel, Marc Lambert and Francis Bach
###   Code supported by Silvere Bonnabel and Marc Lambert  
#############################################################################################

# Project a matrix H=ZZ' Z:d x r on URU'
function [dU,dR,cost]  = proj_lowRank(Z,U,R);
  [d,p]=size(U);
  G=Z'*U;
  dR=G'*G;
  iR=inv(R);
  dU=Z*G*iR-U*dR*iR;
  Q=U'*Z;
  QQ=Z'*Z;
  QG=Z*Q';
  cost=trace(QQ*QQ)-2*trace(QG'*QG)+trace(Q*Q'*Q*Q');
end
