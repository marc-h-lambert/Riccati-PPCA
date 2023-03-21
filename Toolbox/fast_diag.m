############################################################################################
###   Companion code for the article 
###   "Low-rank plus diagonal approximations for Riccati-like matrix differential equations"
###   Authors: Silvere Bonnabel, Marc Lambert and Francis Bach
###   Code supported by Silvere Bonnabel and Marc Lambert  
#############################################################################################

# Fast diagonal of ZZ' when Z: d x r
function [vZ]  = fast_diag(Z);
  [d,r]=size(Z);
  vZ=zeros(d,1);
  for i=1:r;
    vZ=vZ+Z(:,i).^2;
  end
end