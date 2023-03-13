############################################################################################
###   Companion code for the article 
###   "Low-rank plus diagonal approximations for Riccati-like matrix differential equations"
###   Authors: Silvere Bonnabel, Marc Lambert and Francis Bach
###   Code supported by Silvere Bonnabel and Marc Lambert  
#############################################################################################

# Fast computation of the diagonal of XY' where X and Y are of size d x r
function [diag]  = fast_diagXY(X,Y);
  [d,r]=size(X);
  diag=zeros(d,1);
  for i=1:r;
    diag=diag+X(:,i).*Y(:,i);
  end
end