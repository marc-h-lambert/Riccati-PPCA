############################################################################################
###   Companion code for the article 
###   "Low-rank plus diagonal approximations for Riccati-like matrix differential equations"
###   Authors: Silvere Bonnabel, Marc Lambert and Francis Bach
###   Code supported by Silvere Bonnabel and Marc Lambert  
#############################################################################################

function [U,R] = Riccati_lowRank(U,R,dt,A,C,Q,N);
  [d,p]=size(U);
  dU=(eye(d)-U*U')*(A*U+Q*U*inv(R));
  M=U'*A*U*R;
  dR=M+M'+U'*Q*U-R*U'*C'*inv(N)*C*U*R;
  U=retraction_qr(U,dU,dt);
  R=R+dt*(dR);
  R=(R+R')/2;
end