global Ap;
global Bp;
ylim([0 1])
xlim([0 1])

for i=1:size(Ap,1)
  text(Ap(i,2),Ap(i,1), "A_{"+num2str(i)+"}", 'FontSize',18, 'Color','white')
  text(Bp((i),2),Bp((i),1), "B_{"+num2str(i)+"}", 'FontSize',18 ,'Color','Red')
  line(Bp((i),2),Bp((i),1),'Marker', 'O',  'MarkerSize', 14, 'Color', 'white',    'MarkerFaceColor','Red' ) ;
  line(Ap((i),2),Ap((i),1),'Marker', "pentagram"	,  'MarkerSize', 14, 'Color', 'Green',    'MarkerFaceColor','Blue' ) ;

%  line([Ap(i,2), Bp(permut1(i)  ,2)], [Ap(i,1),Bp(permut1 (i),1)],         [0.5,0.5], 'LineStyle', '--','LineWidth',3,'Marker', ">" );
%   line([Ap(i,2), Bp(permut1(i)  ,2)], [Ap(i,1),Bp(permut1 (i),1)],         [0.5,0.5], 'LineStyle', '--','LineWidth',3  );
% 

end
