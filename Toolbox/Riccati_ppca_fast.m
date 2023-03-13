############################################################################################
###   Companion code for the article 
###   "Low-rank plus diagonal approximations for Riccati-like matrix differential equations"
###   Authors: Silvere Bonnabel, Marc Lambert and Francis Bach
###   Code supported by Silvere Bonnabel and Marc Lambert  
#############################################################################################


# U low rank matrix d x p
# R matrix p x p (diagonal)
# s scalar st P=URU'+sI (PPCA decomp)
# sqrtA sqrt of A st sqrtA*sqrtA'=A and sqrtA d x k (k<<d)
# --> if A not symmetric use a sparse structure 
# diagQ diagonal of the process noise matrix d x d 
# C obs matrix m x d
# N obs noise matrix m x m
 
function [U,R,s] = Riccati_ppca_fast(U,R,s,dt,sqrtA,C,diagQ,N);
  [d,p]=size(U);
  
  invN=inv(N);
  K=sqrtA'*U;
  QU=fast_diagU(diagQ,U);
  iRsI=inv(R-s*eye(p));
  dU=sqrtA*K*R+QU-U*(K'*K)*R-U*(U'*QU); # low rank part
  dU=dU+s*sqrtA*K-s*C'*invN*(C*U)*R; # I x s part
  dU=dU-s*U*K'*K+s*U*(U'*(C'*invN*(C*U)*R)); # -UU' x s part
  dU=dU*iRsI;
  U=retraction_qr(U,dU,dt);
  
  M=K'*K*R;
  CU=C*U;
  dR=M+M'+U'*fast_diagU(diagQ,U)-R*CU'*inv(N)*CU*R;
  R=R+dt*(dR);
  R=(R+R')/2;
  
  a1=2*s*sum(fast_diag(sqrtA));
  a2=sum(diagQ);
  #a3=-s*s*trace(invN*C*C');
  a3=-s*s*sum(fast_diagXY(invN*C,C));
  a4=-2*s*sum(fast_diagXY(U*K',sqrtA));
  a5=-sum(fast_diagXY(U,QU));
  a6=s*s*sum(fast_diagXY(invN*CU,CU));
  #a6=s*s*trace(CU'*invN*CU);
  dsigma=0;#(a1+a2+a3+a4+a5+a6)/(d-p);
  s=s+dt*dsigma;
end

##function [U,R,s] = Riccati_ppca(U,R,s,dt,A,C,Q,N);
##  [d,p]=size(U);
##  S=C'*inv(N)*C;
##  dsigma=trace((eye(d)-U*U')*(2*s*A+Q-s*s*S))/(d-p);
##  dU=(eye(d)-U*U')*(A*U*R+Q*U+s*A'*U-s*S*U*R)*inv(R-s*eye(p));
##  dR=U'*A*U*R+R*U'*A*U+U'*Q*U-R*U'*S*U*R;
##  U=retraction_qr(U,dU,dt);
##  R=R+dt*(dR);
##  R=(R+R')/2;
##  s=s+dt*dsigma;
##end
