function [C] = Observation_SwarmDrones(d);
  C=zeros(1,d);
  C(1,d)=1;
  k=randi(d-1);
  b=zeros(1,d);
  b(1,d)=1;
  b(1,k)=-1;
  C=[b; C];
  Na=int8(d/2);
  # for each agents
  for m=1:Na-1;
    # nb neighbors seen =1
    nbAg=randi(1);
    for j=1:nbAg;
      # index of neigbor
      k=randi(Na-1);
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
  b=zeros(1,d);
  b(1,2)=1;
  C=[b; C];
  b=zeros(1,d);
  b(1,1)=1;
  C=[b; C]; %make queen visible
end