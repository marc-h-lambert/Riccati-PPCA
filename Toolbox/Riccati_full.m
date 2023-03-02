############################################################################################
###   Companion code for the article 
###   "Low-rank plus diagonal approximations for Riccati-like matrix differential equations"
###   Authors: Silvere Bonnabel, Marc Lambert and Francis Bach
###   Code supported by Silvere Bonnabel and Marc Lambert  
#############################################################################################

function P = Riccati_full(P,dt,A,C,Q,NN);
  P=P+dt*A*P+dt*P*A'+dt*Q-dt*P*C'*C*P/NN;
  P=(P+P')/2;
end