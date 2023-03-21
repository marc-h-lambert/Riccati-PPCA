############################################################################################
###   Companion code for the article 
###   "Low-rank plus diagonal approximations for Riccati-like matrix differential equations"
###   Authors: Silvere Bonnabel, Marc Lambert and Francis Bach
###   Code supported by Silvere Bonnabel and Marc Lambert  
#############################################################################################

#############################################################################################
# Memory Efficient Riccati with low-rank decomposition P=URU'
# Uncheck all commented code for validation (code commented is not cheap in memory) 
#############################################################################################
# U low rank matrix d x p
# R matrix p x p (diagonal)
# dt integration step 
# sqrtA sqrt of A st sqrtA*sqrtA'=A and sqrtA d x k (k<<d)
# --> if A not symmetric use a sparse structure 
# C obs matrix m x d
# diagQ diagonal of the process noise matrix d x d 
# N obs noise matrix m x m
#############################################################################################
 
function [U,R] = Riccati_lowRank(U,R,dt,sqrtA,C,diagQ,N);
  [d,p]=size(U);
  
  # variables used only for validation
  #Q=diag(diagQ);
  #A=sqrtA*sqrtA';
  #S=C'*inv(N)*C;
  
  # --- compute dU -----
  invR=inv(R);
  K=sqrtA'*U;
  QUiR=fast_diagU(diagQ,U)*invR;
  dU=sqrtA*K+QUiR-U*(K'*K)-U*(U'*QUiR);
  
  # validation : check if dU = dU2
  #dU2=(eye(d)-U*U')*(A*U+Q*U*inv(R))
  
   # --- compute dR -----
  UAUR=K'*K*R;
  CU=C*U;
  dR=UAUR+UAUR'+U'*fast_diagU(diagQ,U)-R*CU'*inv(N)*CU*R;
  
  # validation : check if dR = dR2
  #dR2=U'*A*U*R+R*U'*A'*U+U'*Q*U-R*U'*S*U*R
  
  U=retraction_qr(U,dU,dt);
  R=R+dt*(dR);
  R=(R+R')/2;
end

