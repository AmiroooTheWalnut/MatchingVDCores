function M=Min_Sum_Euclidean_Matching(AX,AY, BX, BY )
nPnts =size(AX,1);
M=1:nPnts;
nV=size(AX,1);
nE=size(AX,1)^2 ; % number of edges
w= 1:nE   ; % zeros(1, nE);

for i=1:nV
  dwAP(i)=NN( AX(i),AY(i), BX,BY  );
  dwBP(i)=NN( BX(i),BY(i), AX,AY );
end

MT1=zeros(nV ,nE); %e_{i,j} is in index 4*(i-1)+(j-1)
MT2=zeros(nV , nE);
k=0;

for i=1:nV
  for j=1:nV

    %
    w( nV*(i-1) +j  )= ((AX(i)-BX(j))^2+ (AY(i)-BY(j))^2  )^0.5;

    MT1( i, nV*(i-1)  + j  ) = 1 ; %edges from A
    MT2( i, nV* (j-1)  +i ) = 1 ;
  end
end

MT3=[MT1;MT2];

f=w;
intcon=1:nE;
b=ones(2*nV ,1 );

lb=zeros(1,nE);
ub=ones(1,nE);


Meq=ones(1,nE);
beq = nV;

x = intlinprog(f,intcon,[] ,[], MT3 , b ,  lb,ub);
%f, x, intcon, b, beq, lb, and ub are vectors, and A and Aeq are matrices.

% min  ( x f ) s.t. 
% x(intcon) are integers
% A⋅x≤b
% Aeq⋅x=beq
% lb≤x≤ub.

for i=1:nV
  for j=1:nV
    if   x( nV*(i-1)+j,1  )  ==1
      M(i)=j;
    end
  end
end
%M=permut1
