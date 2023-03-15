############################################################################################
###   Companion code for the article 
###   "Low-rank plus diagonal approximations for Riccati-like matrix differential equations"
###   Authors: Silvere Bonnabel, Marc Lambert and Francis Bach
###   Code supported by Silvere Bonnabel and Marc Lambert  
#############################################################################################

# Project a matrix H=ZZ' Z:d x r on URU'
function [dU,dR,ds,cost]  = proj_ppca(Z,U,R,s);
  [d,p]=size(U);
  Q=U'*Z;
  QQ=Z'*Z;
  QG=Z*Q';
  ds=(trace(QQ)-trace(U'*Z*Z'*U))/(d-p);
  G=Z'*U;
  dR=G'*G;
  iR=inv(R-s*eye(p));
  dU=Z*G*iR-U*G'*G*iR;
  cost=trace(QQ*QQ)-2*trace(QG'*QG)+trace(Q*Q'*Q*Q');
  cost=cost+ds^2*(d-p)-2*ds*(trace(Z'*Z)-trace(Q*Q'));
end

