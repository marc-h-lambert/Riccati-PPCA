############################################################################################
###   Companion code for the article
###   "Low-rank plus diagonal approximations for Riccati-like matrix differential equations"
###   Authors: Silvere Bonnabel, Marc Lambert and Francis Bach
###   Code supported by Silvere Bonnabel and Marc Lambert
#############################################################################################

addpath Toolbox;
addpath FactorAnalysis;

disp('----------- Test of MEMORY EFFICIENT Riccati flow ------------')
## We integrate the Riccati SDE of the error and the covariance for a simple system


d=1000; # state dim
r=10; # latent dim of P
m=10; # observation dim
ra=10; # latent dim of A

dt=0.01;
Tf=10;
N=Tf/dt;

########## Model of the system

# The dynamic is supposed in the form A=sqrtA*sqrtA' to optimize memory
# --> If A is not symmetric we may use the sparse structure of A to optimize memory
# (using tools of sparse matrix mgt)
sqrtA=zeros(d,ra);
for i=1:ra;
  sqrtA(:,i)=0.1*randn(d,1);
end
  
# Genertate a random observation matrix C of size m x d
C=zeros(m,d);
k=int8(d/10);
for i=1:m;
  Idx=randperm(d,k);
  for j=1:k;
    C(i,Idx)=1;
  endfor
endfor

# The process covariance Q is supposed diagonal
dq=zeros(d,1);
for u=1:d;
  dq(u)=0.1+1*abs(randn);
endfor
NN=0.5*eye(m); %obs noise

########## Compute Pinit of the form URU' (P0 supposed low-rank to be fair with LR decomposition)
TT=zeros(1,N);
O=randn(d,r);
U=orth(O);
UU=U;
R=2*eye(r);
RR=R; 
s=0; # P0 low-rank

%initialisations

%LR KF
U_LR=UU;
R_LR=RR;

%PPCA KF (! watch out !)
U_sI=UU;
R_sI=RR+s*eye(r);

%FA KF (! watch out !)
U_FA=UU;
R_FA=RR;
dpsi_FA=s*ones(d,1);

######## Riccati Low rank #################
fprintf('Low-rank Riccati flow (URU) - %i steps integration  \n',N);
t0=cputime;

for i=1:N-1
  TT(i+1)=dt*i;
  [U_LR,R_LR]=Riccati_lowRank_fast(U_LR,R_LR,dt,sqrtA,C,dq,NN);
end

t=cputime;
timeElapse=(t-t0)*1000;
fprintf('timeElapse = %i ms \n', timeElapse)

######## Riccati PPCA #################
fprintf('Ppca Riccati flow (URU + s(I-UU)) - %i steps integration \n',N);
t0=cputime;

for i=1:N-1
  TT(i+1)=dt*i;
  [U_sI,R_sI,s]=Riccati_ppca_fast(U_sI,R_sI,s,dt,sqrtA,C,dq,NN);
end

t=cputime;
timeElapse=(t-t0)*1000;
fprintf('timeElapse = %i ms \n', timeElapse)

######## Riccati FA #################
fprintf('FA Riccati flow (URU + Psi) - %i steps integration \n',N);
t0=cputime;

for i=1:N-1
  TT(i+1)=dt*i;
  [U_sI,R_sI,s]=Riccati_fa_fast(U_FA,R_FA,dpsi_FA,dt,sqrtA,C,dq,NN);
end

t=cputime;
timeElapse=(t-t0)*1000;
fprintf('timeElapse = %i ms \n', timeElapse)
