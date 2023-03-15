############################################################################################
###   Companion code for the article 
###   "Low-rank plus diagonal approximations for Riccati-like matrix differential equations"
###   Authors: Silvere Bonnabel, Marc Lambert and Francis Bach
###   Code supported by Silvere Bonnabel and Marc Lambert  
#############################################################################################

#############################################################################################
# Memory Efficient Riccati with PPCA decomposition P=URU' + s(I-UU')
# Uncheck all commented code for validation (code commented is not cheap in memory) 
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
 
function [U,R,s] = Riccati_ppca_fast(U,R,s,dt,sqrtA,C,diagQ,N);
  [d,p]=size(U);
  
  # variables used only for validation
  #Q=diag(diagQ);
  #A=sqrtA*sqrtA';
  #S=C'*inv(N)*C;
  
  # --- compute dU -----
  invN=inv(N);
  K=sqrtA'*U;
  QU=fast_diagU(diagQ,U);
  iRsI=pinv(R-s*eye(p));
  dU=sqrtA*K*R+QU-U*(K'*K)*R-U*(U'*QU); # low rank part
  dU=dU+s*sqrtA*K-s*C'*invN*(C*U)*R; # I x s part
  dU=dU-s*U*K'*K+s*U*(U'*(C'*invN*(C*U)*R)); # -UU' x s part
  dU=dU*iRsI;
  
  # validation : check if dU = dU2
  #dU2=(eye(d)-U*U')*(A*U*R+Q*U+s*A'*U-s*S*U*R)*inv(R-s*eye(p))
  
  # --- compute dR -----
  UAUR=K'*K*R;
  CU=C*U;
  dR=UAUR+UAUR'+U'*fast_diagU(diagQ,U)-R*CU'*inv(N)*CU*R;
  
  # validation : check if dR = dR2
  #dR2=U'*A*U*R+R*U'*A'*U+U'*Q*U-R*U'*S*U*R
  
  # --- compute dsigma -----
  a1=2*s*sum(fast_diag(sqrtA));
  a2=sum(diagQ);
  a3=-s*s*sum(fast_diagXY(invN*C,C));
  a4=-2*s*sum(fast_diagXY(U*K',sqrtA));
  a5=-sum(fast_diagXY(U,QU));
  a6=s*s*sum(fast_diagXY(invN*CU,CU));
  dsigma=(a1+a2+a3+a4+a5+a6)/(d-p);
  
  # validation : check if dsigma = dsigma2
  #dsigma2=trace((eye(d)-U*U')*(2*s*A+Q-s*s*S))/(d-p)
  
  # --- project on manifold ----- 
  U=retraction_qr(U,dU,dt);
  R=R+dt*(dR);
  R=(R+R')/2;
  s=s+dt*dsigma;
end
