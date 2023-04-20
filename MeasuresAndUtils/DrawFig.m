function   tmp=DrawFig(M,permut,tit_s)
global Ap;
global Bp;

nexttile
imagesc([0 1],[1 0] ,flip( M(:,:,1) ));
text(0.2, 0.2,  sprintf('%s\n Mean=%5.2f, Max= %5.2f', tit_s ,  mean(M(:,:,1),"all"), max(max(M(:,:,1))) ), ...
    'Color','red', 'BackgroundColor'  ,'blue', 'FontSize',14 );
%title(tit_s)
line([Ap(:,2), Bp( permut(:) ,2)], [Ap(:,1),Bp(permut(:),1)],[0.5,0.5],'LineWidth',2,'Marker', ">", 'color', 'black', ...
    'MarkerSize',18, 'MarkerIndices',[ 2   ], 'LineStyle',"-" ,  'Marker', "*", 'MarkerIndices',[ 1  ]  );
set(gca,'YDir','normal')
make_labels
set(gca,'YDir','normal')



PP=[Ap;Bp] 


%DT=delaunay([Ap(:,2);Bp(:,2) ], [ Ap(:,1);Bp(:,1) ]   )

DT=delaunay(PP(:,2),PP(:,1));

hold on;triplot(DT, [Ap(:,2);Bp(:,2)], [ Ap(:,1);Bp(:,1)],'LineStyle', ':' , 'Color','white');
 


nexttile
imagesc([0 1],[1 0] ,flip( M(:,:,2) ));
%text(0.2, 0.2,  sprintf('%s Mean=%5.2f, Max= %5.2f', tit_s   mean(M(:,:,1),"all"), max(max(M(:,:,1)))), ...
% 'Color','red', 'BackgroundColor'  ,'blue', 'FontSize',16 );
%title(tit_s)
line([Ap(:,2), Bp( permut(:) ,2)], [Ap(:,1),Bp(permut(:),1)],[0.5,0.5],'LineWidth',2,'Marker', ">", 'color', 'black', ...
    'MarkerSize',18, 'MarkerIndices',[ 2   ], 'LineStyle',"-" ,  'Marker', "*", 'MarkerIndices',[ 1  ]  );
set(gca,'YDir','normal')
make_labels
set(gca,'YDir','normal')


end