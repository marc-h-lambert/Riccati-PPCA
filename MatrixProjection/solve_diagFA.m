############################################################################################
###   Companion code for the article 
###   "Low-rank plus diagonal approximations for Riccati-like matrix differential equations"
###   Authors: Silvere Bonnabel, Marc Lambert and Francis Bach
###   Code supported by Silvere Bonnabel and Marc Lambert  
#############################################################################################

# Solve the FA mean square loss ((I-UU') o (I-UU'))^(-1) * d
# where d is supposed to verify Diag (I-UU') M (I-UU')
# (d can be computed using projDiagFA1 or projDiagFA2 
# in function of the structure of M)
function [barphi] = solve_diagFA(diagM,U);
  [d,r]=size(U);
  
  # Compute the matrixes Upsilon and D such that
  # (I-UU') o (I-UU') = I- 2D + Upsilon Upsilon'
  # We have the relation D=Diag(UU') and Upsilon Upsilon'=UU' o UU'
  Upsilon=zeros(d,r*(r+1)/2);
  k=0;
  for j=2:r;
    for i=1:j-1;
      k=k+1;
      Upsilon(:,k)=sqrt(2)*U(:,i).*U(:,j);
    end
  end

  for k=1:r
    Upsilon(:, r*(r-1)/2+k)=U(:,k).*U(:,k);
  end
  
  # For validation check if M1=M2 
  #M1=Upsilon*Upsilon'
  #M2=(U*U').*(U*U')

  diagD=zeros(d,1);
  for i=1:d
    diagD(i)=norm(U(i,:))^2;
  end
  
  % Normal equation solved with memory-cheap Woodbury formulas (if matrixes are invertibles)
    
  if r*(r+1)/2>d ; % we then need a moore penrose pseudo inverse
    % warning moore penrose pseudo inverse may make Kalman diverge
    % warning("r*r > d/2 - fast inverse not used");
    barphi=pinv(eye(d)-2*diag(diagD)+Upsilon*Upsilon')*diagM;
  else
    barpsi=1./(1-2*diagD);
    Z=eye(size(Upsilon)(2))+(Upsilon'.*barpsi')*Upsilon;
    barpsi2=barpsi.*diagM;
    barpsi3=pinv(Z)*Upsilon'*barpsi2;
    barphi=barpsi2-barpsi.*(Upsilon*barpsi3);
  endif
