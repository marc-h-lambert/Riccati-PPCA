############################################################################################
###   Companion code for the article 
###   "Low-rank plus diagonal approximations for Riccati-like matrix differential equations"
###   Authors: Silvere Bonnabel, Marc Lambert and Francis Bach
###   Code supported by Silvere Bonnabel and Marc Lambert  
#############################################################################################

# Fast computation of Diag(v)Z when Z: d x r
function [vZ]  = fast_diagU(v,Z);
  [d,r]=size(Z);
  vZ=zeros(d,r);
  for i=1:r;
    vZ(:,i)=v.*Z(:,i);
  end
end