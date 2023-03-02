############################################################################################
###   Companion code for the article 
###   "Low-rank plus diagonal approximations for Riccati-like matrix differential equations"
###   Authors: Silvere Bonnabel, Marc Lambert and Francis Bach
###   Code supported by Silvere Bonnabel and Marc Lambert  
#############################################################################################

function Y = retraction_qr(X, U, t)
 % It is necessary to call qr_unique rather than simply qr to ensure
       % this is a retraction, to avoid spurious column sign flips.
      if nargin < 3
           Y = qr_unique(X + U);
       else
            Y = qr_unique(X + t*U);
       end
     end
