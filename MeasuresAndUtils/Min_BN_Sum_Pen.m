function M=Min_BN_Sum_Pen(AX,AY, BX, BY )
nPnts =size(AX,1)
M=1:nPnts;
nV=size(AX,1);
nE=size(AX,1)^2 ; % number of edges
w= 1:nE   ; % zeros(1, nE);

BN=eye(nE);

for i=1:nV
    dwAP(i)=NN( AX(i),AY(i), BX,BY  );
    dwBP(i)=NN( BX(i),BY(i), AX,AY );
end

MT1=zeros(nV ,nE); %e_{i,j} is in index 4*(i-1)+(j-1)
MT2=zeros(nV , nE);
k=0;

for i=1:nV
    for j=1:nV
        k=nV*(i-1) +j;
        w( nV*(i-1) +j  )= (((AX(i)-BX(j))^2+ (AY(i)-BY(j))^2  )^0.5)*  (1/dwAP(i) + 1/dwBP(j));


        BN(k,k)=w(k);
        w( nV*(i-1) +j  )= (((AX(i)-BX(j))^2+ (AY(i)-BY(j))^2  )^0.5)*  max(1/dwAP(i) , 1/dwBP(j) );

        MT1( i, nV*(i-1)  + j  ) = 1 ; %edges from A
        MT2( i, nV* (j-1)  +i ) = 1 ;
    end
end



A=[BN,(-1)*ones(nE,1)];
AEQ= [ [MT1;MT2],zeros(2*nV,1)];


f=[zeros(1,nE),1];
intcon=1:nE ;
beq=ones(2*nV ,1 );

lb=zeros(nE+1,1);
ub=[ones(nE,1); 100000];


%Meq=ones(1,nE);
beq = ones(2*nV,1);


x = intlinprog( f, 1:nE,  A,    zeros(nE,1),  ...
    AEQ, beq, lb,ub );
%f, x, intcon, b, beq, lb, and ub are vectors, and A and Aeq are matrices.
% intlinprog(f,intcon,A,b,Aeq,beq,lb,ub,x0,options)
% min  ( x f ) s.t.
% x(intcon) are integers
% A⋅x≤b
% Aeq⋅x=beq
% lb≤x≤ub.
M=zeros(1,nV);
for i=1:nV
    for j=1:nV
        if  abs(  1-  x( nV*(i-1)+j,1  )  ) <0.01 
       
            M(i)=j;
        end
    end
end
%M=permut1
