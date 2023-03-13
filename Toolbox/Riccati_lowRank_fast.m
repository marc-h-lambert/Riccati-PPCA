############################################################################################
###   Companion code for the article 
###   "Low-rank plus diagonal approximations for Riccati-like matrix differential equations"
###   Authors: Silvere Bonnabel, Marc Lambert and Francis Bach
###   Code supported by Silvere Bonnabel and Marc Lambert  
#############################################################################################

# U low rank matrix d x p
# R diag matrix p x p
# sqrtA sqrt of A st sqrtA*sqrtA'=A and sqrtA d x k (k<<d)
# --> if A not symmetric use a sparse structure 
# diagQ diagonal of the process noise matrix d x d 
# C obs matrix m x d
# N obs noise matrix m x m
 
function [U,R] = Riccati_lowRank(U,R,dt,sqrtA,C,diagQ,N);
  [d,p]=size(U);
  invR=inv(R);
  K=sqrtA'*U;
  QUiR=fast_diagU(diagQ,U)*invR;
  dU=sqrtA*K+QUiR-U*(K'*K)-U*(U'*QUiR);
  M=K'*K*R;
  CU=C*U;
  dR=M+M'+U'*fast_diagU(diagQ,U)*invR-R*CU'*inv(N)*CU*R;
  U=retraction_qr(U,dU,dt);
  R=R+dt*(dR);
  R=(R+R')/2;
end

##function [U,R] = Riccati_lowRank(U,R,dt,A,C,phiQ,N);
##  [d,p]=size(U);
##  dU=(eye(d)-U*U')*(A*U+Q*U*inv(R));
##  M=U'*A*U*R;
##  dR=M+M'+U'*Q*U-R*U'*C'*inv(N)*C*U*R;
##  U=retraction_qr(U,dU,dt);
##  R=R+dt*(dR);
##  R=(R+R')/2;
##end