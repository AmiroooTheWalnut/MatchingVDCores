function [M exitflag3] = Min_BN_Pen_Param(AX,AY, BX, BY,...
  r ... %r is the max penalty we are willing to suffer
  )

global exitflag3;
global MaxPenV
global WP;
global Ap ;
global Bp ; 

nPnts =size(AX,1);
M=1:nPnts;


nPnts=size(AX,1);
nE=size(AX,1)^2 ; % number of edges
w= 1:nE   ; % zeros(1, nE);

for i=1:nPnts
  dA(i)=NN( AX(i),AY(i), BX,BY  );
  dB(i)=NN( BX(i),BY(i), AX,AY );
end

MT1=zeros(nPnts ,nE); %e_{i,j} is in index 4*(i-1)+j
MT2=zeros(nPnts , nE);
k=0;

%mdA= min(dA) ;  mdB= min(dB) ;



for i=1:nPnts
  for j=1:nPnts
    k=  nPnts*(i-1)  + j ;% The number of edge (A_i, B_j)
    d=  (((AX(i)-BX(j))^2+ (AY(i)-BY(j))^2  )^0.5) ;
    if d*  max(1/dA(i) , 1/dB(j) ) <=r
      MT1( i, k ) = 1 ; %edges from A
      MT2( j, k ) = 1 ;
    end
  end
end

if (min(max(MT2,[],2))==0 ) || ( min(max(MT1,[],2))==0 )
  exitflag3=-12 % NO matching exists
  return
end

M0=zeros(25,2);
M0(:,1)= (1+mod((1:25),5))/7;
M0(:,2)= (1+round((1:25)/5 ))/7;

[VDV,VDE]=voronoi([AY;BY], [AX;BX])


W4=zeros(2,0)
for i=1:nPnts
  [d0,i0]= min(  ( (AX(i)-BX(:)).^2 + (AY(i)-BY(:)).^2 ).^0.5 )  ;
  W4=[W4;[ 0.3* AY(i)+ 0.7 * BY(i0) ,   0.3* AX(i)+ 0.7 * BX(i0) ]];
  W4=[W4;[ 0.7* AY(i)+ 0.3 * BY(i0) ,   0.7* AX(i)+ 0.3 * BX(i0) ]];

  %   line(W4(i,2),W4(i,1),'Marker', "+"        ,  'MarkerSize', 14, 'Color', 'Green',    'MarkerFaceColor','Blue' ) ;

  [d0,i0]= min(  ( (AX(:)-BX(i)).^2 + (AY(:)-BY(i)).^2 ).^0.5 )  ;
  W4=[W4;[ 0.3* AY(i0)+ 0.7 * BY(i) ,   0.3* AX(i0)+ 0.7 * BX(i) ]];
  W4=[W4;[ 0.7* AY(i0)+ 0.3 * BY(i) ,   0.7* AX(i0)+ 0.3 * BX(i) ]];

end


PP=[Ap;Bp] ;
DT=delaunay([Ap(:,2);Bp(:,2) ], [ Ap(:,1);Bp(:,1) ]   );
W4=zeros(2,0)
figure(9) ; clf; make_labels;hold on
for ind=1:2 %size(DT,1)
  i=DT(ind,1) ; i0=DT(ind,2) ; i1=DT(ind,3);

  W4=[W4;  [0.3* PP(i,:)+ 0.7 * PP(i0,:)] ] ;
  W4=[W4;  [0.7* PP(i,:)+ 0.3 * PP(i0,:)] ] ;
%figure(9);line([PP(i,2), PP(i0,2)] ,  [PP(i,1), PP(i0,1)]  )

  W4=[W4;  [0.3* PP(i0,:)+ 0.7 * PP(i1,:)] ] ;
  W4=[W4;  [0.7* PP(i0,:)+ 0.3 * PP(i1,:)] ] ;
%figure(9);line([PP(i0,2), PP(i1,2)] ,  [PP(i0,1), PP(i1,1)]  )



  W4=[W4;  [0.3* PP(i1,:)+ 0.7 * PP(i,:)] ] ;
  W4=[W4;  [0.7* PP(i1,:)+ 0.3 * PP(i,:)] ] ;
 % figure(9);line([PP(i1,2), PP(i,2)] ,  [PP(i1,1), PP(i,1)]  )

end
WP= W4;
LBmat=zeros(size(WP,1) , nE);

figure(5); clf
make_labels;
hold on ;
text(0.1,0.1, sprintf("r:%4.2f, MaxPen=%4.3f ",r,MaxPenV), ...
  'BackgroundColor'  ,'white', 'color','red', ...
  'FontSize',28 );

for k=1:size(WP,1)
  dA=NN(WP(k,2),WP(k,1), AX,AY);
  dB=NN(WP(k,2),WP(k,1), BX,BY);

  for i=1:nPnts
    for j=1:nPnts

      %line([Ap(i,2), Bp(permut1(i)  ,2)], [Ap(i,1),Bp(permut1 (i),1)],         [0.5,0.5], 'LineStyle', '--','LineWidth',3,'Marker', ">" );
      %   line([Ap(i,2), Bp(permut1(i)  ,2)], [Ap(i,1),Bp(permut1 (i),1)],         [0.5,0.5], 'LineStyle', '--','LineWidth',3  );
      %
      if  (norm( WP(k,:)- Ap(i,:))<=r * dA) && (norm( WP(k,:)- Bp(j,:))<=r * dB )
        line([Ap(i,2), Bp(j  ,2)], [Ap(i,1),Bp(j , 1)],        [0.5,0.5], 'Color','blue',  'LineStyle', '-','LineWidth',1,'Marker', ">" );
        LBmat(k, nPnts*(i-1) +j  )=-1;
      else
        line([Ap(i,2), Bp(j  ,2)], [Ap(i,1),Bp(j , 1)],        [0.5,0.5], 'Color','red', 'LineStyle', '-','LineWidth',1,'Marker', ">" );

      end
    end
  end
  %PossibleEdges((AX,AY, BX, BY, x0,y0); %This function returns 0/1 for the egdes where the egde is admissible for the witness at (x0,yo)
end

b2=zeros(1, size(WP,1));

MT3=[MT1;MT2];

f=ones(1,nE)
intcon=1: nE ;
b=ones(2*nPnts ,1 );

lb=zeros(1,nE);
ub=ones(1,nE);


Meq=ones(1,nE);
beq = nPnts;

[x, fval, exitflag3] = intlinprog(f,intcon,LBmat , b2, MT3 , b ,  lb,ub)
%f, x, intcon, b, beq, lb, and ub are vectors, and A and Aeq are matrices.

% min  ( x f ) s.t.
% x(intcon) are integers
% A⋅x≤b
% Aeq⋅x=beq
% lb≤x≤ub.

if exitflag3~=1
  return
end

for i=1:nPnts
  for j=1:nPnts
    if   x( nPnts*(i-1)+j,1  )  ==1
      M(i)=j;
    end
  end
end

%M=permut
