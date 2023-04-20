


nPnts=8;
Ap=double( zeros(nPnts ,2));
Bp=double( zeros(nPnts ,2));

the=(2*pi/(1 *nPnts) );
Rmat=[[cos(the/2) ,-sin(the/2)];[sin(the/2), cos(the/2)]]
%the=(pi/4 );

for i=1:nPnts ;
  Ap(i,1)=sin(  1* (i-1) *the ) ; %The Y value
  Ap(i,2)=cos( 1* (i-1) *the ) ;%The X value
  %   Bp(i,1)=sin( 2* (i-1) *the +the ) ;  %The Y value
  %   Bp(i,2)=cos( 2* (i-1) *the +the ) ;
end
%Ap= 0.5 + [  0.4 *Ap;  0.3*Ap*Rmat]
Bp=   [0.7 *Ap*Rmat]

Ap=0.5 +0.5*Ap;
Bp=0.5 +0.5*Bp;

%Bp=rescale(Bp)
%
% for i=1:nPnts
%   Ap(i,2)=(0.5+i-1)  *1/nPnts  ;
%   Bp(i,2)=(0.5+ i-1) *1/nPnts ;
%   Ap(i,1)=0.1   ;
%   Bp(i,1)=0.1
% end
% Ap(4,2)=Ap(1,2); Ap(4,1)=1-Ap(1,1);
% Bp(4,2)=Bp(1,2)+0.3; Bp(4,1)=Ap(4,1);

%  Ap=rand(nPnts,2)
%  Bp=rand(nPnts,2)
%
%    Ap(:,:)=[[0.5 0.1]; [0.45 0.3]; [0.7 0.3]];
%    Bp(:,:)=[[0.5 0.9]; [0.65 0.3] ;[0.71 0.35] ];





NP=false;
if(NP==true)
  nPnts=11;


  Bp=[ [0.5,0.86603];  [2.071,0]; [1.3355,2.31315] ]
  Ap=[ [1,0] ; [1.0355,1.79354];[2.671,0]]
  t=2*pi/3
  R=[[(cos(t)), -1*sin( t) ; (sin( t)),cos(t)]]
  Ap=[Ap;Ap*R;Ap*R*R]/6+0.5
  Bp=[Bp;Bp*R; Bp*R*R]/6+0.5

  Ap=[Ap;[0.69,0.9]; [0.45,0.9] ]
  Bp=[Bp;[0.55,0.9];  [0.35,0.9]]
end
%
%
% Ap=[[0.5 0.1]; [0.5, 0.3] ; [0.9 0.8]]
% Bp=[[0.55 0.6];[0.55, 0.7]; [0.9 0.9]]
% minX1=min([Ap;Bp](:,2))
%  m2=max([Ap;Bp],[],'all')
%Ap= rescale(Ap,m1,m2);Bp= rescale(Bp,m1,m2);

  Ap=rand(nPnts,2)
  Bp=rand(nPnts,2)
% figure
% line([Ap(:,2), Bp( : ,2)], [Ap(:,1),Bp(:,1)],[0.5,0.5],'LineWidth',5,'Marker', ">", 'color', 'black', ...
%   'MarkerSize',18, 'MarkerIndices',[ 2   ], 'LineStyle',"-" ,  'Marker', "*", 'MarkerIndices',[ 1  ]  );

