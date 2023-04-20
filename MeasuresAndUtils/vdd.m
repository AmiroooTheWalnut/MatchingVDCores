q17=zeros(n,n,2);
 
global imdA 
global imdB


 
for i=1:n
  for j=1:n
    % dd=zeros(nPnts,2);
    x=j/n ; y=i/n;
    %   [x,y]=ginput(1)
    

    dA=  ((x-Ap(:,2)).^2          + (y-Ap(:,1)).^2 ).^0.5  ;

    dB=  ((x-Bp(:,2)).^2 + (y-Bp(:,1)).^2 ).^0.5;
    [mdA, imdA]=min(dA) ; [mdB ,imdB] = min(dB);

    if (mdA<=mdB)
      %   q17(i,j,:)= [mdA imdA];
      pen= max(  norm( [y,x]- Bp(   permut_Min_BN_Pen_Param(imdA),:  ) ) /mdB , ...
        norm( [y,x]- Ap(imdA,:) ) / mdA);
      q17(i,j,:)=[ pen, imdA] ;
      %    sprintf("case 1, ")

    else %Nearest point is in B
      i0= find(  permut_Min_BN_Pen_Param==imdB );
      pen= max(  norm( [y,x]- Bp(   imdB, :  ) ) /mdB , ...
        norm( [y,x]- Ap(i0,: ) ) /mdA);
      %   sprintf("case 2, %d %d ", imdB, i0)

      q17(i,j,:)=[ pen, i0 ] ;% % imdB+nPnts] ;

    end
  end

end


figure(7); clf
% 
% 
% DrawFig( q17   , permut_Min_BN_Pen_Param , ...
%   sprintf(  "Any matching with Max Ben < %4.3f", r) );
% V_Pen(Ap(:,2),Ap(:,1),Bp(:,2), Bp(:,1), permut_Min_BN_Pen_Param)
%  
% voronoi( [Ap(:,2);Bp(:,2)], [Ap(:,1);Bp(:,1)]  );



