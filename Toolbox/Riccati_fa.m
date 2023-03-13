############################################################################################
###   Companion code for the article 
###   "Low-rank plus diagonal approximations for Riccati-like matrix differential equations"
###   Authors: Silvere Bonnabel, Marc Lambert and Francis Bach
###   Code supported by Silvere Bonnabel and Marc Lambert  
#############################################################################################

function [U,R,Psi] = Riccati_fa(U,R,Psi,dt,A,C,Q,N);
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
end