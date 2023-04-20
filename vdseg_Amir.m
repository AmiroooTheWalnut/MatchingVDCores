
global n ;
global nPnts  ;

global MaxPenV %Max Pen on a vertex 
global Ap;
global Bp;

global permut_Min_Sum_Euc ;  %This is the matching that minimizes  $$\sum_{p_i,q_j\in M}  distance(p_i,q_j) $$
global permut_Min_Sum_Pen ;  %This is the matching that minimizes  $$\sum_{p_i,q_j\in M}  w(s_i,g_j)$,  where $(s_i, q_j)=\max( s_i/ \delta(s_i, G), g_j/delta( g_j, S)$
%global permut_Min_Sum_Squares ;
global permut_Min_BN_Euc
global permut_Min_BN_Pen;
global permut_Min_BN_Pen_Param;
global permut_Min_BN_Sum_Pen;

global q_Min_BN_Sum_Pen;

global q_Added_Pen ;  %The cost is the round-trip walk to gym and store, and the VD is w resapect to the segment (matched pair)  that takes mi distance to reach.
%   The matching is the same one used in permut_Min_Sum_Pen

global q_Min_Sum_Squares;

global q_Min_Sum_Euc ;
global q_Min_Sum_Pen ;
global q_Min_BN_Pen ;
global q_Min_BN_Euc ; 
global q_Min_BN_Pen_Param;
global exitflag3;
global WP;

addpath("MeasuresAndUtils\")


n=190;

rng(8);

% CreateInput;

numShops=6;
numGyms=6;
shops=rand(numShops,2);
gyms=rand(numGyms,2);
Ap=shops;
Bp=gyms;


nPnts =size(Ap,1)

[permut_Min_BN_Pen, MaxPenV]= Min_BN_Pen( Ap(:,2),Ap(:,1),Bp(:,2), Bp(:,1) );
sprintf("computed Bottleneck Pen   \n\n",'color','red')


for r=max(1,MaxPenV):0.1:(2*MaxPenV+1) 
  [permut_Min_BN_Pen_Param,exitflag] =Min_BN_Pen_Param( Ap(:,2),Ap(:,1),Bp(:,2), Bp(:,1), r );
  if exitflag==1 
    break ;
  end 
end


 
permut_Min_Sum_Euc=     Min_Sum_Euclidean_Matching( Ap(:,2),Ap(:,1),Bp(:,2), Bp(:,1) );
disp("computed Min Sum Euc \n\n")

permut_Min_Sum_Pen=     Min_Sum_Penalty_Matching( Ap(:,2),Ap(:,1),Bp(:,2), Bp(:,1) );
disp("computed Min Sum Euc \n\n")


permut_Min_BN_Euc= Min_BN_Euc( Ap(:,2),Ap(:,1),Bp(:,2), Bp(:,1) ); %Bottleneck
sprintf("computed Bottleneck Euc   \n\n")

[permut_Min_BN_Pen, MaxPenV]= Min_BN_Pen( Ap(:,2),Ap(:,1),Bp(:,2), Bp(:,1) );
sprintf("computed Bottleneck Pen   \n\n",'color','red')


permut_Min_BN_Sum_Pen = Min_BN_Sum_Pen( Ap(:,2),Ap(:,1),Bp(:,2), Bp(:,1) );
sprintf("computed BN  Sum Pen Pen   \n\n",'color','red')
%The weight of an edge is the sum of vertices's weights




q_Min_Sum_Euc=      zeros(n,n,2) ;
q_Min_BN_Pen=      zeros(n,n,2) ;
q_Min_Sum_Pen=     zeros(n,n,2) ;
q_Min_BN_Euc=         zeros(n,n,2) ;
q_Min_BN_Param=         zeros(n,n,2) ;

sprintf("computed Bottleneck Penalties   \n\n")

q_Added_Pen  =        zeros(n,n,2) ;

%
%
%   PP_Min_Sum_EucA =  zeros(1,nPnts);
%  PP_Min_Sum_PenB ;   zeros(1,nPnts);
%  PP_Min_Sum_SquaresB=zeros(1,nPnts);
%  PP_Min_BN_Euc=zeros(1,nPnts);
%  PP_Min_BN_Pen=zeros(1,nPnts);




%figure(1) ; clf(1)
%setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);

for i=1:n
    for j=1:n
        % dd=zeros(nPnts,2);
        x=j/n ; y=i/n;

        dA=  ((x-Ap(:,2)).^2          + (y-Ap(:,1)).^2 ).^0.5  ;





        dB=  ((x-Bp(:,2)).^2 + (y-Bp(:,1)).^2 ).^0.5;
        [mdA imdA]=min(dA) ; [mdB imdB] = min(dB);
        [u v]=min( max(    [dA /mdA, dB/mdB ] , [], 2)) ;
        q0(i,j,[1 2])=[u v];



        dB=  ((x-Bp(permut_Min_Sum_Euc(:),2)).^2 + (y-Bp(permut_Min_Sum_Euc(:),1)).^2 ).^0.5;
        mdA=min(dA) ; mdB = min(dB);
        [u v]=min( max(    [dA /mdA, dB/mdB ] , [], 2)) ;
        q_Min_Sum_Euc(i,j,[1 2])=[u v];



        dB=  ((x-Bp(permut_Min_Sum_Pen(:),2)).^2 + (y-Bp(permut_Min_Sum_Pen(:),1)).^2 ).^0.5;
        [u v]=min( max(    [dA/mdA , dB/mdB ] , [], 2)) ;
        q_Min_Sum_Pen(i,j,[1 2])=[u v];



        dB=  ((x-Bp(permut_Min_BN_Sum_Pen(:),2)).^2 + (y-Bp( permut_Min_BN_Sum_Pen(:),1)).^2 ).^0.5;
        [u v]=min((dA+dB)/(mdA+mdB));  %Q ADDED _ NOTE THE DIFFERENCE
        q_Added_Pen(i,j,[1 2])=[u v];



        dB=  ((x-Bp(permut_Min_BN_Pen(:),2)).^2 + (y-Bp(permut_Min_BN_Pen(:),1)).^2 ).^0.5;
        [u v]=min( max(    [dA/mdA , dB/mdB ] , [], 2)) ;
        q_Min_BN_Pen(i,j,[1 2])=[u v];



        dB=  ((x-Bp(permut_Min_BN_Pen_Param(:),2)).^2 + (y-Bp(permut_Min_BN_Pen_Param(:),1)).^2 ).^0.5;
        [u v]=min( max(    [dA/mdA , dB/mdB ] , [], 2)) ;
        q_Min_BN_Pen_Param(i,j,[1 2])=[u v];


    end


end

f=figure(1) ; clf(1)
%f = figure(WindowKeyPressFcn=@figureCallback);

setappdata(gcf, 'SubplotDefaultAxesLocation', [0, 0, 1, 1]);
t = tiledlayout(3,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';
t.TileSpacing = 'none';


%-------
%DrawFig(q0, 1:nPnts,"Arbitrary matching")


DrawFig(q_Min_Sum_Euc, permut_Min_Sum_Euc , "Min Sum Euclidean Matching)")
V_Pen(Ap(:,2),Ap(:,1),Bp(:,2), Bp(:,1), permut_Min_Sum_Euc)

DrawFig(q_Min_Sum_Pen ,permut_Min_Sum_Pen,"Min Sum of Penelties  " );
V_Pen(Ap(:,2),Ap(:,1),Bp(:,2), Bp(:,1), permut_Min_Sum_Pen)


DrawFig(q_Min_BN_Pen  , permut_Min_BN_Pen,  "Min Bottlenck of Penelty " );
V_Pen(Ap(:,2),Ap(:,1),Bp(:,2), Bp(:,1), permut_Min_BN_Pen)

% % DrawFig( q_Added_Pen  , permut_Min_Sum_Pen , ...
% %     "Minimize round trip  gym+show normalized by round to nearests" );
% 
% V_Pen(Ap(:,2),Ap(:,1),Bp(:,2), Bp(:,1), permut_Min_BN_Sum_Pen)


 

q2=q_Min_BN_Pen_Param;  ;
for i=1:n
  for j=1:n
    q2(i,j,1)=max(q2(i,j,1),2);
  end
end
 
DrawFig( q2   , permut_Min_BN_Pen_Param , ...
    sprintf(  "Any matching with Max Ben < %4.3f", r) );
V_Pen(Ap(:,2),Ap(:,1),Bp(:,2), Bp(:,1), permut_Min_BN_Pen_Param)



vdd; figure(1)
DrawFig( q17   , permut_Min_BN_Pen_Param , ...
    sprintf(  "Matching to nearest site < %4.3f", r) );
V_Pen(Ap(:,2),Ap(:,1),Bp(:,2), Bp(:,1), permut_Min_BN_Pen_Param)


figure(8)


DrawFig( q17   , permut_Min_BN_Pen_Param , ...
    sprintf(  "Matching to nearest site < %4.3f", r) );
V_Pen(Ap(:,2),Ap(:,1),Bp(:,2), Bp(:,1), permut_Min_BN_Pen_Param)

figure(9) 


datacursormode on
% make_labels;
% colormap;
%
%

%------
%
%
% MaxPenalty=max(max(q(:,:,1)));
%
% [xx,yy] = meshgrid(0:1/n:1-1/n , 0:1/n:1-1/n  );
%
% % linkaxes([ad1, ad2,ad3], 'xy')
%
% figure(3)
% alpha 0.5
%
%
% surf(xx,yy,q(:,:,2)) ; hold on; mesh(xx,yy,q2(:,:,2));
%
% colorbar
% figure(4)
% alpha 0.5
%
% mesh(xx,yy,q(:,:,1),q(:,:,2))
% set(gca,'YDir','normal')
%
%
%
%
%
% %linkaxes([ae1, ae2], 'xyz')
%
% %figure ; imshow( imresize( RGB, 7) );
% s="VDexample"+date+".mat";
% save(s,'A','B','f3', '-v7.3','-nocompression')

%imwrite(RGB, "VDexample"+date+".jpg");
%savefig(date+".fig");

writematrix(Bp,"Bp"+date+".txt");

