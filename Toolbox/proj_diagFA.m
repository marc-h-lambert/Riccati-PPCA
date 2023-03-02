############################################################################################
###   Companion code for the article 
###   "Low-rank plus diagonal approximations for Riccati-like matrix differential equations"
###   Authors: Silvere Bonnabel, Marc Lambert and Francis Bach
###   Code supported by Silvere Bonnabel and Marc Lambert  
#############################################################################################

function [res,cost] = proj_diagFA(G,U);
  d=length(G(:,1));
  r=length(U(1,:));
  %we implement barh=diag(G*G');
  barh=zeros(d,1);
  for i=1:d
    barh(i)=norm(G(i,:))^2;
  end

  %we implement baruh=diag(U*(U'*G)*G');
  UG=(U'*G)*G';
  baruh=zeros(d,1);
  for i=1:d
    baruh(i)=U(i,:)*UG(:,i);
  end

  %we now compute Lambda


  lambda=zeros(d,1);
  Q=U'*G;
  for i=1:d
    %lambda(i)=trace(U'*G*G'*U*U(i,:)'*U(i,:));
    B=U(i,:)*Q;
    lambda(i)=trace(B'*B);
  end


  %Computation of matrix Upsilon
    Upsilon=zeros(d,r*(r+1)/2);
    k=0;
  for j=2:r;
    for i=1:j-1;
      k=k+1;
      Upsilon(:,k)=sqrt(2)*U(:,i).*U(:,j);
    end
  end

  for k=1:r
    Upsilon(:, r*(r-1)/2+k)=U(:,k).*U(:,k);
  end

  V=zeros(d,1);
  for i=1:d
    V(i)=norm(U(i,:))^2;
  end
  D=diag(V);

  %big matrix with an inversion
  if r^2/2>d % we then need a moore penrose pseudo inverse
    barphi=pinv(eye(d)-2*D+Upsilon*Upsilon')*(barh-2*baruh+lambda);
    res=barphi;
  else
    %we need to code vectorially the following
    Psi=diag(1./diag(eye(d)-2*D));
    %Inv=Psi-Psi*Upsilon*inv(eye(size(Upsilon)(2))+Upsilon'*Psi*Upsilon)*Upsilon'*Psi;
    %barphi=Inv*(barh-2*baruh+lambda)  ;
    %en vectoriel
    barpsi=1./diag(eye(d)-2*D);
    barpsi2=barpsi.*(barh-2*baruh+lambda) ;
    barpsi3=pinv(eye(size(Upsilon)(2))+Upsilon'*Psi*Upsilon)*Upsilon'*barpsi2;
    barphi=barpsi2-barpsi.*(Upsilon*barpsi3);
    res=barphi;
  endif

  cout00=norm(barh-barphi)^2+4*baruh'*barphi-2*barphi'*D*barphi-2*lambda'*barphi;
  cout00=cout00+ barphi'*Upsilon*Upsilon'*barphi;
  cost=cout00-norm(barh)^2;
end
