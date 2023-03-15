############################################################################################
###   Companion code for the article 
###   "Low-rank plus diagonal approximations for Riccati-like matrix differential equations"
###   Authors: Silvere Bonnabel, Marc Lambert and Francis Bach
###   Code supported by Silvere Bonnabel and Marc Lambert  
#############################################################################################

#############################################################################################
# Memory Efficient Riccati with FA decomposition P=URU' + Psi
# Uncheck all commented code for validation (code commented is not cheap in memory) 
#############################################################################################
# U low rank matrix d x p
# R matrix p x p (diagonal)
# diagPsi vector of dim d correponding to the diagonal of Psi
# dt integration step 
# sqrtA sqrt of A st sqrtA*sqrtA'=A and sqrtA d x k (k<<d)
# --> if A not symmetric use a sparse structure 
# C obs matrix m x d
# diagQ diagonal of the process noise matrix d x d 
# N obs noise matrix m x m
#############################################################################################
 

function [U,R,diagPsi] = Riccati_fa_fast(U,R,diagPsi,dt,sqrtA,C,diagQ,N);  
  # define intermediate variables
  [d,p]=size(U);  
  invN=inv(N);
  invR=inv(R);
  K=sqrtA'*U;
  UAU=K'*K;
  QU=fast_diagU(diagQ,U);
  CU=C*U;
  PsiU=fast_diagU(diagPsi,U);
  PsiC=fast_diagU(diagPsi,C');
  PsiA=fast_diagU(diagPsi,sqrtA);
  CPsiU=C*PsiU;
  
##  # VALIDATION variables used only for validation
##  Q=diag(diagQ);
##  A=sqrtA*sqrtA';
##  Psi=diag(diagPsi);
##  I=eye(d);
##  Piuu=I-U*U';
##  S=C'*inv(N)*C;
##  M1=A*Psi+Psi*A'+Q-Psi*S*Psi;
  
  # --- compute dPsi -----
  d=fast_diagFA1(sqrtA,PsiA,U);
  d=d+fast_diagFA1(PsiA,sqrtA,U);
  d=d+fast_diagFA2(diagQ,U);
  d=d-fast_diagFA1(PsiC*invN,PsiC,U);
  ddiagPsi=solve_diagFA(d,U);
  
##  # VALIDATION : check if d = mhm and ddiagPsi = ddiagPsi2
##  iPUU2=pinv(Piuu.*Piuu);
##  mhm=diag(Piuu*(M1)*Piuu); # must be equal to d
##  ddiagPsi2=iPUU2*mhm; # must be equal to ddiagPsi
    
  # --- compute dU -----
  DU=fast_diagU(diagQ-ddiagPsi,U);
  UUPsiSU=U*(CPsiU'*invN*CU);
  PsiSU=PsiC*(invN*CU);
  sAPsiUiR=(sqrtA'*PsiU)*invR;
  PsiAUiR=PsiA*(sqrtA'*U)*invR;
  UUAPsiUiR=U*K'*sAPsiUiR;
  PsiSPsiUiR=PsiC*invN*CPsiU*invR;
  
  dU=sqrtA*sAPsiUiR+PsiAUiR+QU*invR;
  dU=dU-fast_diagU(ddiagPsi,U)*invR-PsiSPsiUiR;
  dU=dU-UUAPsiUiR-U*(U'*PsiAUiR);
  dU=dU-U*(U'*DU*invR)+U*(U'*PsiSPsiUiR);
  dU=dU+sqrtA*K-PsiSU-U*UAU+UUPsiSU;
  
##  # VALIDATION : check if dU = dU2
##  dU2=Piuu*(M1-diag(ddiagPsi))*U*inv(R)+Piuu*(A*U-Psi*S*U);
  
  # --- compute dR -----
  RUSPsiU=R*CU'*inv(N)*(C*PsiU);
  UAPsiU=K'*(sqrtA'*PsiU);
  
  dR=UAU*R+R'*UAU'-R*CU'*inv(N)*CU*R; 
  dR=dR-RUSPsiU-RUSPsiU';
  dR=dR+UAPsiU+UAPsiU';
  dR=dR-CPsiU'*invN*CPsiU; # terms between U' and U
  dR=dR+U'*DU; # terms between U' and U
  
##  # VALIDATION : check if dR = dR2
##  M2=R*U'*S*Psi*U;
##  M3=U'*A*U*R;
##  dR2=U'*(M1-diag(ddiagPsi))*U-M2-M2'-R*U'*S*U*R+M3+M3';
    
  # --- project on manifold ----- 
  U=retraction_qr(U,dU,dt);
  R=R+dt*(dR);
  R=(R+R')/2;
  diagPsi=diagPsi+dt*ddiagPsi;
end
