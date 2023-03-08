############################################################################################
###   Companion code for the article
###   "Low-rank plus diagonal approximations for Riccati-like matrix differential equations"
###   Authors: Silvere Bonnabel, Marc Lambert and Francis Bach
###   Code supported by Silvere Bonnabel and Marc Lambert
#############################################################################################

addpath Toolbox;

'----------- Diagonal + low rank : Test of Riccati flow ------------'
## We integrate the Riccati SDE of the error and the covariance for a simple system
# with A=0 and P0=URU'+s(I-UU')

# Can change d or r or Cvariable or SQ
# The run is senstive to how SQ and U are generared
# Cvariable=0: 1D problem with fixed C
# Cvariable=1: 2D problem with variable C

d=150;
r=10;

%Simulation du système

dt=0.01;
Tf=10;
N=Tf/dt;

% Tout le KF peut être fait hors ligne (à moins qu'on modifie le graphe en fonction des distances).

########## Model of the system
Cvariable=0;
A=zeros(d,d);

#if Cvariable ==0;
  # Observation matrix the same at each t
  C=zeros(d,d);
  for i=1:d-1;
    C(i,i+1)=1;
  end
  C=-eye(d)+C;
  C(d,d)=1;
  C;
#endif

dq=zeros(d,1);
for u=1:d;
  dq(u)=0.1+1*abs(rand);
endfor
dq=dq;
#dq=randn(d,1);
#SQ=diag(1:d)/d;
#SQ=0.1*diag(rand(d));
Q=10*0.2*diag(dq);
#SQ=eye(d);
#Q=0.1*SQ*SQ'; %process noise
SQ=chol(Q);
dq=zeros(d,1);
for u=1:d;
  dq(u)=0.1+1*abs(randn);
endfor
dq=dq;
NN=2*diag(dq);

Q=100*0.01*diag(dq);
#SQ=eye(d);
#Q=0.1*SQ*SQ'; %process noise
SQ=chol(Q);
NN=10*0.2; %obs noise





########## Compute Pinit
TT=zeros(1,N);
O=randn(d,r);
UU=orth(O);
U=UU;%;eye(d)(:,1:r); % the goal is to generate a stiefel random matrix to avoid alignment with the diagonal which would favor our methods
RR=2*eye(r); %initial parameters for all filters


R=RR;
s=1*0.1;

P_init=U*R*U' + s*eye(d);  %+s*(eye(d)-0*U*U'); % this is the common initial P

#initial error
eX=randn(d);

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


%FA KF v2
U_FA2=UU;
R_FA2=RR;
Psi_FA2=s*eye(d);
P_FA2=P_init;

eX_KF=eX;
eX_LR=eX;
eX_sI=eX;
eX_FA=eX;

########## Compute Kalman
tt=zeros(4,N);
ee=zeros(4,N);
for i=1:N-1
  TT(i+1)=dt*i;

# Change Matrix C
if mod(i*dt,5)==0 && Cvariable!=0;
  '---regenerate Matrix C -----'
  i
  C=zeros(1,d);
  C(1,d)=1;
  k=randi(d-1)
  b=zeros(1,d);
  b(1,d)=1;
  b(1,k)=-1;
  C=[b; C];
  Na=int8(d/2);
  # for each agents
  for m=1:Na-1;
    # nb neighbors between 1 and 3
    nbAg=randi(3);
    for j=1:nbAg;
      # index of neigbor
      k=randi(Na);
      if k!=m;
        b=zeros(1,d);
        b(1,2*(m-1)+1)=-1;
        b(1,2*(k-1)+1)=1;
        C=[b; C];
        b=zeros(1,d);
        b(1,2*(m-1)+2)=-1;
        b(1,2*(k-1)+2)=1;
        C=[b; C];
      endif
    end
  end
  C;
endif



% Full KF
S=C'*inv(NN)*C;
P_KF=Riccati_full(P_KF,dt,A,C,Q,NN);
dX_KF=(A-P_KF*S)*eX_KF;
eX_KF=eX_KF+dX_KF*dt;

% PPCA KF
[U_sI,R_sI,s]=Riccati_ppca(U_sI,R_sI,s,dt,A,C,Q,NN);
P_sI=U_sI*R_sI*U_sI'+s*(eye(d)-U_sI*U_sI');
dX_sI=(A-P_sI*S)*eX_sI;
eX_sI=eX_sI+dX_sI*dt;

% Low_rank KF
[U_LR,R_LR]=Riccati_lowRank(U_LR,R_LR,dt,A,C,Q,NN);
P_LR=U_LR*R_LR*U_LR';
dX_LR=(A-P_LR*S)*eX_LR;
eX_LR=eX_LR+dX_LR*dt;

##% FA KF v2 (with vectorize form)
##[U_FA2,R_FA2,Psi_FA2]=Riccati_fa(U_FA2,R_FA2,Psi_FA2,dt,A,C,Q,NN);
##P_FA2=U_FA2*R_FA2*U_FA2'+Psi_FA2;
##dX_FA=(A-P_FA2*S)*eX_FA;
##eX_FA=eX_FA+dX_FA*dt;
##norm(P_KF-P_FA2)/norm(P_KF); # (the two versions must give the same error)

# Validation of FA KF with successive projection
G1=SQ;
phi1=diag(proj_diagFA(G1,U_FA));
U1= (eye(d)-U_FA*U_FA')*(G1*G1'-phi1)*U_FA*inv(R_FA);
D1=U_FA'*(G1*G1'-phi1)*U_FA;

G2=(U_FA*R_FA*U_FA'+psi_FA)*C'/sqrt(NN);
phi2=diag(proj_diagFAneg(G2,U_FA));
U2=(eye(d)-U_FA*U_FA')*(-G2*G2'-phi2)*U_FA*inv(R_FA);
D2=U_FA'*(-G2*G2'-phi2)*U_FA;

psi_FA=psi_FA+dt*(phi1+phi2);
dU=(U1+U2);
U_FA=retraction_qr(U_FA,dU,dt);
R_FA=R_FA+dt*(D1+D2);
R_FA=(R_FA+R_FA')/2;
P_FA=U_FA*R_FA*U_FA'+psi_FA;
dX_FA=(A-P_FA*S)*eX_FA;
eX_FA=eX_FA+dX_FA*dt;
norm(P_KF-P_FA)/norm(P_KF); #(the two versions must give the same error)

tt(1,i+1)= norm(P_KF-P_LR)/norm(P_KF);
tt(2,i+1)= norm(P_KF-P_sI)/norm(P_KF);
tt(3,i+1)= norm(P_KF-P_FA)/norm(P_KF);

#ee(1,i+1)= norm(eX_LR);
#ee(2,i+1)= norm(eX_sI);
#ee(3,i+1)= norm(eX_FA);
#ee(4,i+1)= norm(eX_KF);

ee(1,i+1)= norm(eX_KF-eX_LR);
ee(2,i+1)= norm(eX_KF-eX_sI);
ee(3,i+1)= norm(eX_KF-eX_FA);
ee(4,i+1)= norm(eX_KF-eX_KF);



end

'psi_FA'
det(psi_FA)

figure(1),clf
linewidth=2;
fontsize=30;
plot(TT(1,:),tt(1,:),'-.',"linewidth",linewidth,TT(1,:),tt(2,:),"linewidth",linewidth,TT(1,:),tt(3,:),"linewidth",linewidth)
axis([0 Tf 0 1]);
legend('Low-rank','PPCA','FA','Location','northwest');
xlabel ("time (s)");
ylabel ("error on P");
title('||P-proj(P)|| / ||P||','fontsize',20)
%h=get(gcf, "currentaxes");
%set(h, "fontsize", fontsize, "linewidth", linewidth);
print "-S200,200" -dpdf -color cov_XP3.pdf

figure(2)
plot(TT(1,:),ee(1,:),'-.',"linewidth",linewidth,TT(1,:),ee(2,:),"linewidth",linewidth,TT(1,:),ee(3,:),"linewidth",linewidth,TT(1,:),ee(4,:),"linewidth",linewidth)
axis([0 Tf 0 max(max(ee))+2]);
xlabel ("time (s)");
ylabel ("error on X");
legend('Low-rank','PPCA','FA','KF','Location','southwest')
%h=get(gcf, "currentaxes");
%set(h, "fontsize", fontsize, "linewidth", linewidth);
title('||dX||')
print "-S200,200" -dpdf -color err_XP3.pdf
