############################################################################################
###   Companion code for the article 
###   "Low-rank plus diagonal approximations for Riccati-like matrix differential equations"
###   Authors: Silvere Bonnabel, Marc Lambert and Francis Bach
###   Code supported by Silvere Bonnabel and Marc Lambert  
#############################################################################################

# Helper method for factor analysis decomposition
# Compute Diag of (I-UU')M(I-UU') when M=XY' where X and Y are of dim d x p (p<<<d)
function [diagPhi] = fast_diagFA1(X,Y,U);
  [d,r]=size(U);
  h=fast_diagXY(X,Y);
  hu1=fast_diagXY(U,Y*(X'*U));
  Z=U*(U'*Y);
  hu2=fast_diagXY(X,Z);
  Lambda=fast_diagXY(U*(U'*X),Z);
  diagPhi=h-hu1-hu2+Lambda;
end
