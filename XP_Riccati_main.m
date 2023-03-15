############################################################################################
###   Companion code for the article
###   "Low-rank plus diagonal approximations for Riccati-like matrix differential equations"
###   Authors: Silvere Bonnabel, Marc Lambert and Francis Bach
###   Code supported by Silvere Bonnabel and Marc Lambert
#############################################################################################

addpath Toolbox;
addpath FactorAnalysis;

disp('----------- Diagonal + low rank : Test of Riccati flow ------------')
## We integrate the Riccati SDE of the error and the covariance for a simple system
# with A=0 and P0=URU'+s(I-UU')

disp('--------- XP with dim 200 and latent dimension 8 ------------')
d=200;
r=8;

dt=0.01;
Tf=10;
N=Tf/dt;

########## Model of the system
ra=1;
sqrtA=zeros(d,ra);
A=sqrtA*sqrtA';
C=Observation_SwarmDrones(d);

dq=zeros(d,1);
for u=1:d;
  dq(u)=0.1+1*abs(randn);
endfor
Q=diag(dq);
SQ=diag(sqrt(dq));
NN=2; %obs noise

########## Compute Pinit
TT=zeros(1,N);
O=randn(d,r);
UU=orth(O);
U=UU;%;eye(d)(:,1:r); % the goal is to generate a stiefel random matrix to avoid alignment with the diagonal which would favor our methods
RR=2*eye(r); %initial parameters for all filters
R=RR;
s=0;

P_init=U*R*U' + s*eye(d);  %+s*(eye(d)-0*U*U'); % this is the common initial P


########## Initialisations

% KF
P_KF=P_init;

%LR KF
U_LR=UU;
R_LR=RR;
P_LR=P_init;

%PPCA KF (! watch out !)
U_sI=UU;
R_sI=RR+s*eye(r);
P_sI=P_init;

%FA KF
U_FA=UU;
R_FA=RR;
dpsi_FA=diag(s*eye(d));
P_FA=P_init;

#initial error
eX=randn(d);
eXX=eX;
eX_KF=eX;
eX_LR=eX;
eX_sI=eX;
eX_FA=eX;

########## Compute Kalman
tt=zeros(3,N);
ee=zeros(4,N);
for i=1:N-1
  TT(i+1)=dt*i;
  
  # Change Matrix C to test a changing graph
  if mod(i*dt,12)==0;
    disp('---regenerate Matrix C -----')
    C=Obsrervation_SwarmDrones(d);
  endif
  
  if mod(i*dt,1)==0;
    fprintf('step-%i \n',i);
  endif
  
  % Full KF
  if mod(i*dt,1)==0;
    disp('step Full KF')
  endif
  S=C'*inv(NN)*C;
  P_KF=Riccati_full(P_KF,dt,sqrtA,C,dq,NN);
  dX_KF=(A-P_KF*S)*eX_KF;
  eX_KF=eX_KF+dX_KF*dt;
  
  % Low_rank KF
  if mod(i*dt,1)==0;
    disp('step Riccati-LR')
  endif
  [U_LR,R_LR]=Riccati_lowRank_fast(U_LR,R_LR,dt,sqrtA,C,dq,NN);
  P_LR=U_LR*R_LR*U_LR';
  dX_LR=(A-P_LR*S)*eX_LR;
  eX_LR=eX_LR+dX_LR*dt;
  
  % PPCA KF
  if mod(i*dt,1)==0;
    disp('step Riccati-PPCA')
  endif
  [U_sI,R_sI,s]=Riccati_ppca_fast(U_sI,R_sI,s,dt,sqrtA,C,dq,NN);
  P_sI=U_sI*R_sI*U_sI'+s*(eye(d)-U_sI*U_sI');
  dX_sI=(A-P_sI*S)*eX_sI;
  eX_sI=eX_sI+dX_sI*dt;
  
  % FA KF 
  if mod(i*dt,1)==0;
    disp('step Riccati-FA')
  endif
  [U_FA,R_FA,dpsi_FA]=Riccati_fa_fast(U_FA,R_FA,dpsi_FA,dt,sqrtA,C,dq,NN);
  P_FA=U_FA*R_FA*U_FA'+diag(dpsi_FA);
  dX_FA=(A-P_FA*S)*eX_FA;
  eX_FA=eX_FA+dX_FA*dt;
  
  tt(1,i+1)= norm(P_KF-P_LR)/norm(P_KF);
  tt(2,i+1)= norm(P_KF-P_sI)/norm(P_KF);
  tt(3,i+1)= norm(P_KF-P_FA)/norm(P_KF);
  
  ee(1,i+1)= norm(eX_KF-eX_LR);
  ee(2,i+1)= norm(eX_KF-eX_sI);
  ee(3,i+1)= norm(eX_KF-eX_FA);
  ee(4,i+1)= norm(eX_KF-eX_KF);

endfor

figure(1),clf
linewidth=2;
fontsize=30;
plot(TT(1,:),tt(1,:),'-.',"linewidth",linewidth,TT(1,:),tt(2,:),"linewidth",linewidth,TT(1,:),tt(3,:),"linewidth",linewidth)
axis([0 Tf 0 1]);
legend('Low-rank','PPCA','FA','Location','north');
xlabel ("time (s)");
ylabel ("error on P");
title('||P-P_{KF}|| / ||P_{KF}||','fontsize',20)
h=get(gcf, "currentaxes");
set(h, "fontsize", fontsize, "linewidth", linewidth);
print "-S200,200" -dpdf -color cov_XP2.pdf

figure(2)
plot(TT(1,:),ee(1,:),'-.',"linewidth",linewidth,TT(1,:),ee(2,:),"linewidth",linewidth,TT(1,:),ee(3,:),"linewidth",linewidth)
axis([0 Tf 0 max(max(ee))+2]);
xlabel ("time (s)");
ylabel ("error on X w.r.t. Full KF");
legend('Low-rank','PPCA','FA','Location','northwest')
h=get(gcf, "currentaxes");
set(h, "fontsize", fontsize, "linewidth", linewidth);
title('||X-X_{KF}||')
print "-S200,200" -dpdf -color err_XP2.pdf


'--------- XP with dim 200 and latent dimension 50 ------------'
yl=max(max(ee))+2;
r=50;

########## Compute Pinit
TT=zeros(1,N);
O=randn(d,r);
UU=orth(O);
U=UU;%;eye(d)(:,1:r); % the goal is to generate a stiefel random matrix to avoid alignment with the diagonal which would favor our methods
RR=2*eye(r); %initial parameters for all filters
R=RR;
s=0;

P_init=U*R*U' + s*eye(d);  %+s*(eye(d)-0*U*U'); % this is the common initial P


%initialisations
% KF
P_KF=P_init;

%PPCA KF (! watch out !)
U_sI=UU;
R_sI=RR+s*eye(r);
P_sI=P_init;

%LR KF
U_LR=UU;
R_LR=RR;
P_LR=P_init;

%FA KF
U_FA=UU;
R_FA=RR;
psi_FA=s*eye(d);
P_FA=P_init;

#initial error
eX=eXX;
eX_KF=eX;
eX_LR=eX;
eX_sI=eX;
eX_FA=eX;

########## Compute Kalman
tt=zeros(3,N);
ee=zeros(4,N);
for i=1:N-1
  TT(i+1)=dt*i;
  
  # Change Matrix C to test a changing graph
  if mod(i*dt,12)==0;
    disp('---regenerate Matrix C -----')
    C=Obsrervation_SwarmDrones(d);
  endif
  
  if mod(i*dt,1)==0;
    fprintf('step-%i \n',i);
  endif
  
  % Full KF
  if mod(i*dt,1)==0;
    disp('step Full KF')
  endif
  S=C'*inv(NN)*C;
  P_KF=Riccati_full(P_KF,dt,sqrtA,C,dq,NN);
  dX_KF=(A-P_KF*S)*eX_KF;
  eX_KF=eX_KF+dX_KF*dt;
  
  % Low_rank KF
  if mod(i*dt,1)==0;
    disp('step Riccati-LR')
  endif
  [U_LR,R_LR]=Riccati_lowRank_fast(U_LR,R_LR,dt,sqrtA,C,dq,NN);
  P_LR=U_LR*R_LR*U_LR';
  dX_LR=(A-P_LR*S)*eX_LR;
  eX_LR=eX_LR+dX_LR*dt;
  
  % PPCA KF
  if mod(i*dt,1)==0;
    disp('step Riccati-PPCA')
  endif
  [U_sI,R_sI,s]=Riccati_ppca_fast(U_sI,R_sI,s,dt,sqrtA,C,dq,NN);
  P_sI=U_sI*R_sI*U_sI'+s*(eye(d)-U_sI*U_sI');
  dX_sI=(A-P_sI*S)*eX_sI;
  eX_sI=eX_sI+dX_sI*dt;
  
  % FA KF 
  if mod(i*dt,1)==0;
    disp('step Riccati-FA')
  endif
  [U_FA,R_FA,dpsi_FA]=Riccati_fa_fast(U_FA,R_FA,dpsi_FA,dt,sqrtA,C,dq,NN);
  P_FA=U_FA*R_FA*U_FA'+diag(dpsi_FA);
  dX_FA=(A-P_FA*S)*eX_FA;
  eX_FA=eX_FA+dX_FA*dt;
  
  tt(1,i+1)= norm(P_KF-P_LR)/norm(P_KF);
  tt(2,i+1)= norm(P_KF-P_sI)/norm(P_KF);
  tt(3,i+1)= norm(P_KF-P_FA)/norm(P_KF);
  
  ee(1,i+1)= norm(eX_KF-eX_LR);
  ee(2,i+1)= norm(eX_KF-eX_sI);
  ee(3,i+1)= norm(eX_KF-eX_FA);
  ee(4,i+1)= norm(eX_KF-eX_KF);

endfor

figure(3),clf
linewidth=2;
fontsize=30;
plot(TT(1,:),tt(1,:),'-.',"linewidth",linewidth,TT(1,:),tt(2,:),"linewidth",linewidth,TT(1,:),tt(3,:),"linewidth",linewidth)
axis([0 Tf 0 1]);
legend('Low-rank','PPCA','FA','Location','north');
xlabel ("time (s)");
ylabel ("error on P");
title('||P-P_{KF}|| / ||P_{KF}||','fontsize',20)
h=get(gcf, "currentaxes");
set(h, "fontsize", fontsize, "linewidth", linewidth);
print "-S200,200" -dpdf -color cov_XP1.pdf

figure(4)
plot(TT(1,:),ee(1,:),'-.',"linewidth",linewidth,TT(1,:),ee(2,:),"linewidth",linewidth,TT(1,:),ee(3,:),"linewidth",linewidth)
axis([0 Tf 0 yl]);
xlabel ("time (s)");
ylabel ("error on X w.r.t. Full KF");
legend('Low-rank','PPCA','FA','Location','northwest')
h=get(gcf, "currentaxes");
set(h, "fontsize", fontsize, "linewidth", linewidth);
title('||X-X_{KF}||')
print "-S200,200" -dpdf -color err_XP1.pdf
