############################################################################################
###   Companion code for the article 
###   "Low-rank plus diagonal approximations for Riccati-like matrix differential equations"
###   Authors: Silvere Bonnabel, Marc Lambert and Francis Bach
###   Code supported by Silvere Bonnabel and Marc Lambert  
#############################################################################################

currentDir=fileparts(mfilename('fullpath'));
addpath(fullfile(currentDir,'Toolbox'));
addpath(fullfile(currentDir,'MatrixProjection'));

disp('----------- Diagonal + low rank : Test of Projection ------------');
## We project a matrix H on the hyperplane of the LR/FA/PPCA manifold
# The hyperplane is evaluated around a point U0 R0 Psi0

d=1000000; # dimension (rank of diagonal Psi)
r=10; # rank of low-rank part URU'
p=100; # rank of the target matrix H

######## Hypothesis on H ###########
G=randn(d,p);
%G=G.*G;
%H=G*G' #-> G is the factor of matrix H we want to project
%det(H)

######## Hypothesis on local hyperplane ###########
U0=eye(d)(:,1:r); # a stiefel matrix corresponding to the first r canonical vectors   
R0=2*diag(1:r); # a diagonal low rank matrix

### ML initialize  H0=GG'+Diag(0) such that all methods start from the same local point ###
### If Diag term is not null, not clear how to estimate the local point for LR and PPCA ??
barpsi0=zeros(d,1);
s0=sum(barpsi0)/d;

######## Project LR #################
disp('Projection on -- low rank -- manifold around U0 R0');
t0=cputime;
[U,R,cost]  = proj_lowRank(G,U0,R0);
t=cputime;
timeElapse=(t-t0)*1000;
fprintf('timeElapse with "low-rank" method = %i ms \n', timeElapse)
#error=sqrt(cost);#1-cost/trace(QQ*QQ);
#fprintf('cost low-rank (normalized) = %i\n', error)

##
######## Project PPCA #################
disp('Projection on -- ppca -- manifold around U0 R0 s0');
t0=cputime;
[U,R,s,cost]  = proj_ppca(G,U0,R0,s0);
t=cputime;
timeElapse=(t-t0)*1000;
fprintf('timeElapse with "ppca" method = %i ms \n', timeElapse)
#error=sqrt(cost)#1-cost/trace(QQ*QQ);
#fprintf('cost ppca (normalized) = %i\n', error)

######## Project FA #################
disp('Projection on -- fa -- manifold around U0 R0 Psi0');
t0=cputime;
[U,R,barpsi]=proj_fa(G,U0,R0,barpsi0);
t=cputime;
timeElapse=(t-t0)*1000;
fprintf('timeElapse with "fa" method = %i ms \n', timeElapse);
#error=sqrt(cost)#1-cost/trace(QQ*QQ);
#fprintf('cost fa (normalized) = %i\n', error)