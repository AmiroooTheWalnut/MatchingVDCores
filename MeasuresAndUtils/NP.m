
close all; figure 
 

 


 Bp=[ [0.5,0.86603];  [2.071,0]; [1.3355,2.31315] ]
 Ap=[ [1,0] ; [1.0355,1.79354];[2.671,0]]
 t=2*pi/3 
 R=[[(cos(t)), -1*sin( t) ; (sin( t)),cos(t)]]
Ap=[Ap;Ap*R;Ap*R*R]
Bp=[Bp;Bp*R; Bp*R*R]

line([Ap(:,2), Bp( : ,2)], [Ap(:,1),Bp(:,1)],[0.5,0.5],'LineWidth',5,'Marker', ">", 'color', 'black', ...
'MarkerSize',18, 'MarkerIndices',[ 2   ], 'LineStyle',"-" ,  'Marker', "*", 'MarkerIndices',[ 1  ]  );
