############################################################################################
###   Companion code for the article 
###   "Low-rank plus diagonal approximations for Riccati-like matrix differential equations"
###   Authors: Silvere Bonnabel, Marc Lambert and Francis Bach
###   Code supported by Silvere Bonnabel and Marc Lambert  
#############################################################################################

function [U,R,s] = Riccati_ppca(U,R,s,dt,A,C,Q,N);
  [d,p]=size(U);
  dsigma=trace((eye(d)-U*U')*(2*s*A+Q-s*s*C'*inv(N)*C))/(d-p);
  dU=(eye(d)-U*U')*(A*U*R+Q*U+s*A'*U-s*C'*inv(N)*C*U*R)*inv(R-s*eye(p));
  dR=U'*A*U*R+R*U'*A*U+U'*Q*U-R*U'*C'*inv(N)*C*U*R;
  U=retraction_qr(U,dU,dt);
  R=R+dt*(dR);
  R=(R+R')/2;
  s=s+dt*dsigma;
end