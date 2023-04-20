function M=V_Pen(AX,AY, BX, BY ,permut )
nPnts =size(AX,1);
M=1:nPnts;
nV=size(AX,1);
global WP;

PenVA=zeros(1,nV);
PenVB=zeros(1,nV);

for i=1:nV
  dwAP(i)=NN( AX(i),AY(i), BX,BY  );
  dwBP(permut(i))=NN( BX(permut(i)),BY(permut(i)), AX,AY );
  d=(  ( AX(i)-BX(permut(i)) )^2+    (AY(i)-BY(permut(i)))^2  )^0.5;
  PenVA(i)= d/ dwAP(i);
  PenVB(permut(i))=d/ dwBP(permut(i)) ;
  text(AX(i),AY(i), sprintf(" %4.2f", PenVA(i)) , ...
    'BackgroundColor'  ,'white',  'Color','blue',   'FontSize',10 );
  text(BX(permut(i)),BY(permut(i)), sprintf(" %4.2f", PenVB(permut(i))) , ...
    'BackgroundColor'  ,'white', 'color','red', ...
    'FontSize',11 );


end


[a i ]=max(PenVA) ;
text(AX(i),AY(i), sprintf("   %4.2f",  PenVA(i)) ,  'Color','white', 'BackgroundColor'  ,'blue', 'FontSize',15 );

[a i ]=max(PenVB) ;
text(BX(i),BY(i), sprintf("  %4.2f", PenVB(i)) ,  'Color','black', 'BackgroundColor'  ,'red', 'FontSize',15 );

%
 DT=delaunay([AX;BX], [AY;BY] )
 hold on;triplot(DT,  [AX;BX], [AY;BY]  ,'k-','Color','white', 'LineStyle', ':');



for i=1:size(WP,1)
  text(WP(i,2),WP(i,1), sprintf("* %3d",i) , ...
    'BackgroundColor'  ,'white',  'Color','blue',   'FontSize',13 );
end
end

