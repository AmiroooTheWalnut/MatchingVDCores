function M=Min_Sum_Penalty_Matching(AX,AY, BX, BY,...
    r ... %r is the max penalty we are willing to suffer
    )

global exitflag3;

nPnts =size(AX,1)
M=1:nPnts


figure; make_labels;
nV=size(AX,1);
nE=size(AX,1)^2 ; % number of edges
w= 1:nE   ; % zeros(1, nE);

for i=1:nV
    dA(i)=NN( AX(i),AY(i), BX,BY  );
    dB(i)=NN( BX(i),BY(i), AX,AY );
end

MT1=zeros(nV ,nE); %e_{i,j} is in index 4*(i-1)+(j-1)
MT2=zeros(nV , nE);
k=0;

mdA= min(dA) ;  mdB= min(dB) ;


for i=1:nV
    for j=1:nV
        w( nV*(i-1) +j  )= (((AX(i)-BX(j))^2+ (AY(i)-BY(j))^2  )^0.5)*  max(1/dA(i) , 1/dB(j) );
        if w( nV*(i-1) +j  )> r
            w( nV*(i-1) +j  )=10000
        end
        MT1( i, nV*(i-1)  + j  ) = 1 ; %edges from A
        MT2( i, nV* (j-1)  +i ) = 1 ;
    end
end



M0=zeros(25,2);
M0(:,1)= (1+mod((1:25),5))/7;
M0(:,2)= (1+round((1:25)/5 ))/7;


WP= M0;
LBmat=zeros(size(WP,1) , nE);

figure
make_labels
for k=1:size(WP,1)

    dA=NN(WP(k,2),WP(k,1), AX,AY);
    dB=NN(WP(k,2),WP(k,1), BX,BY);


    for i=1:nV
        for j=1:nV
hold on
            %line([Ap(i,2), Bp(permut1(i)  ,2)], [Ap(i,1),Bp(permut1 (i),1)],         [0.5,0.5], 'LineStyle', '--','LineWidth',3,'Marker', ">" );
            %   line([Ap(i,2), Bp(permut1(i)  ,2)], [Ap(i,1),Bp(permut1 (i),1)],         [0.5,0.5], 'LineStyle', '--','LineWidth',3  );
            %
            if  (norm( WP(k,:)- Ap(i,:))<=r * dA) && (norm( WP(k,:)- Bp(j,:))<=r * dB )
                line([Ap(i,2), Bp(j  ,2)], [Ap(i,1),Bp(j , 1)],        [0.5,0.5], 'Color','blue',  'LineStyle', '-','LineWidth',1,'Marker', ">" );
                 LBmat(k, nV*(i-1) +j  )=-1;

            else
               %                 line([Ap(i,2), Bp(j  ,2)], [Ap(i,1),Bp(j , 1)],        [0.5,0.5], 'Color','red', 'LineStyle', '-','LineWidth',1,'Marker', ">" );

            end
        end
    end
    %PossibleEdges((AX,AY, BX, BY, x0,y0); %This function returns 0/1 for the egdes where the egde is admissible for the witness at (x0,yo)
end

b2=zeros(1, size(WP,1))

MT3=[MT1;MT2];

f=w
intcon=1: nE ;
b=ones(2*nV ,1 );

lb=zeros(1,nE);
ub=ones(1,nE);


Meq=ones(1,nE);
beq = nV;

[x, fval, exitflag3] = intlinprog(f,intcon,LBmat , b2, MT3 , b ,  lb,ub)
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

%M=permut
