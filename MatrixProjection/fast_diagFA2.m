############################################################################################
###   Companion code for the article 
###   "Low-rank plus diagonal approximations for Riccati-like matrix differential equations"
###   Authors: Silvere Bonnabel, Marc Lambert and Francis Bach
###   Code supported by Silvere Bonnabel and Marc Lambert  
#############################################################################################

# Helper method for factor analysis decomposition
# Compute Diag of (I-UU') M (I-UU') when M=Diag(diagPsi)
function [diagPhi] = fast_diagFA2(diagPsi,U);
  [d,r]=size(U);
  h=diagPsi;
  PsiU=fast_diagU(diagPsi,U);
  hu=fast_diagXY(PsiU,U);
  UPsiU=U'*PsiU;
  Lambda=fast_diagXY(U,U*UPsiU');
  diagPhi=h-2*hu+Lambda;
end
