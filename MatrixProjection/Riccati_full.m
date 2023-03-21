############################################################################################
###   Companion code for the article 
###   "Low-rank plus diagonal approximations for Riccati-like matrix differential equations"
###   Authors: Silvere Bonnabel, Marc Lambert and Francis Bach
###   Code supported by Silvere Bonnabel and Marc Lambert  
#############################################################################################

#############################################################################################
# Standard Riccati equations 
#############################################################################################
# P matrix d x d 
# dt integration step 
# sqrtA sqrt of A st sqrtA*sqrtA'=A and sqrtA d x k (k<<d)
# --> if A not symmetric use a sparse structure 
# C obs matrix m x d
# diagQ diagonal of the process noise matrix d x d 
# N obs noise matrix m x m
#############################################################################################
 
function P = Riccati_full(P,dt,sqrtA,C,diagQ,NN);
  CP=C*P;
  P=P+2*dt*sqrtA*(sqrtA'*P)+dt*diag(diagQ)-dt*CP'*inv(NN)*CP;
  P=(P+P')*0.5;
end