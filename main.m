clc
clear
%rng(2157)
rng(8)
%rng(7)%interesting, one circle in the middle

numShops=3;
numGyms=3;

global shops;
global gyms;
global isOnLine

isOnLine=false;

%\/\/\/ random?
shops=rand(numShops,2);
gyms=rand(numGyms,2);
if isOnLine==true
    for i=1:numShops
        shops(i,2)=rand(1,1)*0.001;
    end
    for i=1:numGyms
        gyms(i,2)=rand(1,1)*0.001;
    end
end
%^^^ random?
%shops=[0.906,0.55;0.98,0.22;0.88,0.01];
%gyms=[0.51,0.12;0.35,0.58;0.65,0.32];
shops=[0.906000000000000,0.550000000000000;0.980000000000000,0.220000000000000;0.880000000000000,0.0100000000000000];
gyms=[0.488532110091743,0.223241590214067;0.383027522935780,0.538226299694190;0.687355783340648,0.426575066892618];

global shopsPOIs;
global gymsPOIs;

shopsPOIs=cell(1,numShops);
gymsPOIs=cell(1,numShops);
% Lroi=cell(1,n);
for i=1:numShops
    shopsPOIs{1,i} = images.roi.Point(gca,'Position',shops(i,:),'Color','b');
    addlistener(shopsPOIs{1,i},'ROIMoved',@allevents);
end
for i=1:numGyms
    gymsPOIs{1,i} = images.roi.Point(gca,'Position',gyms(i,:),'Color','g');
    addlistener(gymsPOIs{1,i},'ROIMoved',@allevents);
end
run(shops,gyms);

function allevents(src,evt)
global shopsPOIs;
global gymsPOIs;
global isOnLine
n=size(shopsPOIs,2);
shopPositions(n,2)=0;
for i=1:n
    shopPositions(i,:)=shopsPOIs{1,i}.Position;
end
n=size(gymsPOIs,2);
gymPositions(n,2)=0;
for i=1:n
    gymPositions(i,:)=gymsPOIs{1,i}.Position;
end
if isOnLine==true
    for i=1:size(shopPositions,1)
        shopPositions(i,2)=rand(1,1)*0.001;
    end
    for i=1:size(gymPositions,1)
        gymPositions(i,2)=rand(1,1)*0.001;
    end
end
figure(1);
evname = evt.EventName;
axesHandlesToChildObjects = findobj(gca, 'Type', 'line');
if ~isempty(axesHandlesToChildObjects)
    delete(axesHandlesToChildObjects);
end
switch(evname)
    case{'MovingROI'}
        run(shopPositions,gymPositions);
    case{'ROIMoved'}
        run(shopPositions,gymPositions);
end
end

function run(shops,gyms)
global shopsPOIs;
global gymsPOIs;
global isOnLine
figure(1)
clf
hold on
if isOnLine==true
    ylim([-0.5,0.5])
    xlim([-1,2])
else
    ylim([0,1])
    xlim([0,1])
end
daspect([1,1,1])
% \/\/\/ GET CORES AND REGIONS
drawForbiddenAreas(shops,gyms);
[SVCells,GVCells]=calcCores(shops,gyms);%DONE
drawCore(SVCells,GVCells,shops,gyms);%DONE
% ^^^ GET CORES AND REGIONS
numShops=size(shops,1);
for i=1:numShops
    shopsPOIs{1,i} = images.roi.Point(gca,'Position',shops(i,:),'Color','b');
    %addlistener(shopsPOIs{1,i},'MovingROI',@allevents);
    addlistener(shopsPOIs{1,i},'ROIMoved',@allevents);
    text(shops(i,1)+0.01,shops(i,2)+0.01,strcat('Shop',num2str(i)))
end
numGyms=size(gyms,1);
for i=1:numGyms
    gymsPOIs{1,i} = images.roi.Point(gca,'Position',gyms(i,:),'Color','g');
    %addlistener(shopsPOIs{1,i},'MovingROI',@allevents);
    addlistener(gymsPOIs{1,i},'ROIMoved',@allevents);
    text(gyms(i,1)+0.01,gyms(i,2)+0.01,strcat('Gym',num2str(i)))
end

if isOnLine==true
    % SAMPLE FROM SPACE
    points=linspace(-1,2,1000);
    [musts,forbiddens]=assignForbiddensMusts1D(points,shops,gyms);
else
    x=-2:0.05:2;
    y=-2:0.05:2;
    [X,Y] = meshgrid(x,y);
    points(2,size(x,2)*size(y,2))=0;
    counter=1;
    for i=1:size(X,1)
        for j=1:size(X,2)
            points(:,counter)=[X(i,j),Y(i,j)];
            counter=counter+1;
        end
    end
    [musts,forbiddens]=assignForbiddensMusts2D(points,shops,gyms);
end
distances(size(shops,1),size(gyms,1))=0;
for i=1:numShops
    for j=1:numGyms
        distances(i,j)=pdist([shops(i,:);gyms(j,:)]);
    end
end
[f,intcon,A,b,Aeq,beq,lb,ub]=createILP(musts,forbiddens,distances);
x = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub)
if isempty(x)==0
    x=decodeSolution(x,numShops,numGyms,true);
    sum(sum(x.*distances))
end
disp('!!!')
end

function x=decodeSolution(x,numShops,numGyms,isReport)
x=reshape(x,[numShops,numGyms])';
if isReport==true
    for i=1:numShops
        for j=1:numGyms
            if x(i,j)==1
                fprintf('Shop %i is connected to Gym %i\n', i,j)
            end
        end
    end
end
end

function [f,intcon,A,b,Aeq,beq,lb,ub]=createILP(musts,forbiddens,distances)
%f=ones(size(musts,1)*size(musts,2),1);
f=reshape(distances',[size(musts,1)*size(musts,2),1]);
intcon=1:size(musts,1)*size(musts,2);
A(2*size(musts,1)*size(musts,2),size(musts,1)*size(musts,2))=0;
b(2*size(musts,1)*size(musts,2),1)=0;
counter=1;
for i=1:size(musts,1)
    for j=1:size(musts,2)
        col=((i-1)*size(musts,2))+j;
        A(counter,col)=-1;
        b(counter,1)=-musts(i,j);
        counter=counter+1;
    end
end
for i=1:size(musts,1)
    for j=1:size(musts,2)
        col=((i-1)*size(musts,2))+j;
        A(counter,col)=1;
        b(counter,1)=1-forbiddens(i,j);
        counter=counter+1;
    end
end
counter=1;
Aeq(size(musts,1)+size(musts,2),size(musts,1)*size(musts,2))=0;
beq(size(musts,1)+size(musts,2),1)=0;
for i=1:size(musts,1)
    for j=1:size(musts,2)
        col=((i-1)*size(musts,2))+j;
        Aeq(counter,col)=1;
    end
    beq(counter,1)=1;
    counter=counter+1;
end
for j=1:size(musts,2)
    for i=1:size(musts,1)
        col=((i-1)*size(musts,2))+j;
        Aeq(counter,col)=1;
    end
    beq(counter,1)=1;
    counter=counter+1;
end
lb(size(musts,1)*size(musts,2),1)=0;
ub(size(musts,1)*size(musts,2),1)=0;
ub(:,1)=1;
end

function [musts,forbiddens]=assignForbiddensMusts2D(points,shops,gyms)
forbiddens(size(shops,1),size(gyms,1))=0;
musts(size(shops,1),size(gyms,1))=0;
for i=1:size(points,2)
    %points(1,i)=0.3;
    isInShopCores=isInCore(points(:,i),shops);
    isInGymCores=isInCore(points(:,i),gyms);
    isInShopForbiddens=isInForbidden(points(:,i),shops);
    isInGymForbiddens=isInForbidden(points(:,i),gyms);
    for m=1:size(shops,1)
        for n=1:size(gyms,1)
            if isInShopCores(1,m)==1 && isInGymCores(1,n)==1
                musts(m,n)=1;
            end
        end
    end
    for m=1:size(shops,1)
        for n=1:size(gyms,1)
            if isInShopForbiddens(1,m)==1 && isInGymCores(1,n)==1
                forbiddens(m,n)=1;
                 %if forbiddens(2,1)==1
                 %    disp('!!!')
                 %end
            end
        end
    end
    for n=1:size(shops,1)
        for m=1:size(gyms,1)
            if isInGymForbiddens(1,m)==1 && isInShopCores(1,n)==1
                forbiddens(n,m)=1;
%                 if forbiddens(1,3)==1
%                     disp('!!!')
%                 end
            end
        end
    end
%     scatter(points(1,i),0,100)
    
%     if sum(sum(forbiddens))>5
%         disp('!!!')
%     end
    %disp('!!!')
end
%disp('!!!')
end

function [musts,forbiddens]=assignForbiddensMusts1D(points,shops,gyms)
forbiddens(size(shops,1),size(gyms,1))=0;
musts(size(shops,1),size(gyms,1))=0;
for i=1:size(points,2)
    %points(1,i)=0.3;
    isInShopCores=isInCore([points(1,i);0],shops);
    isInGymCores=isInCore([points(1,i);0],gyms);
    isInShopForbiddens=isInForbidden([points(1,i);0],shops);
    isInGymForbiddens=isInForbidden([points(1,i);0],gyms);
    for m=1:size(shops,1)
        for n=1:size(gyms,1)
            if isInShopCores(1,m)==1 && isInGymCores(1,n)==1
                musts(m,n)=1;
            end
        end
    end
    for m=1:size(shops,1)
        for n=1:size(gyms,1)
            if isInShopForbiddens(1,m)==1 && isInGymCores(1,n)==1
                forbiddens(m,n)=1;
%                 if forbiddens(1,3)==1
%                     disp('!!!')
%                 end
            end
        end
    end
    for n=1:size(shops,1)
        for m=1:size(gyms,1)
            if isInGymForbiddens(1,m)==1 && isInShopCores(1,n)==1
                forbiddens(n,m)=1;
%                 if forbiddens(1,3)==1
%                     disp('!!!')
%                 end
            end
        end
    end
%     scatter(points(1,i),0,100)
    
%     if sum(sum(forbiddens))>5
%         disp('!!!')
%     end
    %disp('!!!')
end
%disp('!!!')
end

function isInCores=isInCore(point,POIs)
isInCores(1,size(POIs,1))=0;
for i=1:size(POIs,1)
    isInCoreTemp=1;
    for j=1:size(POIs,1)
        if i~=j
            if pdist([point';POIs(i,:)])>pdist([point';POIs(j,:)])/2
                isInCoreTemp=0;
                break;
            end
        end
    end
    isInCores(1,i)=isInCoreTemp;
    if isInCoreTemp==1
        %scatter(POIs(i,1),POIs(i,2),200);
        %scatter(point(1,1),0,200,'filled');
        %disp('!!!')
    end
end

%disp('!!!')
end

function isInForbiddens=isInForbidden(point,POIs)
isInForbiddens(1,size(POIs,1))=0;
for i=1:size(POIs,1)
    for j=1:size(POIs,1)
        if i~=j
            if pdist([point';POIs(j,:)])<pdist([point';POIs(i,:)])/2
                isInForbiddens(1,i)=1;
                break;
            end
        end
    end
    if isInForbiddens(1,i)==1
        %scatter(POIs(i,1),POIs(i,2),200);
        %scatter(point(1,1),0,200,'filled');
        %disp('!!!')
    end
end
%disp('!!!')
end

function drawForbiddenAreas(shops,gyms)
h=voronoi(shops(:,1),shops(:,2));
set(h, 'Color', 'b')

h=voronoi(gyms(:,1),gyms(:,2));
set(h, 'Color', 'g')

%legend('Shop POI','Shop forbidden','Gym POI','Gym forbidden')

dt = delaunayTriangulation(shops);
[SSV,SSR] = voronoiDiagram(dt);
dt = delaunayTriangulation(gyms);
[GSV,GSR] = voronoiDiagram(dt);

SVCells=calcNeighbors(SSV,SSR,shops,false);
GVCells=calcNeighbors(GSV,GSR,gyms,false);

[allRegions,allCircles,allIntersections]=getAllRegions(SVCells,GVCells,SSR,GSR);

%disp('!')
end

function [allRegions,allCircles,allIntersections]=getAllRegions(SVCells,GVCells,SSR,GSR)
allRegions=cell(1,4);
allRegions{1,1}='Region name';
allRegions{1,2}='Region owners';
allRegions{1,3}='Region forbidden to owners';
allRegions{1,4}='Arcs';
allCircles=cell(1,6);
allCircles{1,1}='Circle index';
allCircles{1,2}='Circle owner';
allCircles{1,3}='Circle forbidden to owner';
allCircles{1,4}='Circle center';
allCircles{1,5}='Circle radius';
allCircles{1,6}='Is excluded';
allIntersections=cell(1,7);
allIntersections{1,1}='Circle 1 index';
allIntersections{1,2}='Circle 2 index';
allIntersections{1,3}='Point';
allIntersections{1,4}='Is checked';
allIntersections{1,5}='Angle to circle 1';
allIntersections{1,6}='Angle to circle 2';
allIntersections{1,7}='Angle to circle candidate';
counter=1;
for i=1:size(SVCells,1)-1
    for j=1:size(SVCells{i+1,2},2)
        allCircles{counter+1,1}=counter;
        allCircles{counter+1,2}=strcat('S',num2str(SVCells{i+1,2}(1,j)));
        allCircles{counter+1,3}=strcat('S',num2str(i));
        allCircles{counter+1,4}=SVCells{i+1,3}(:,j);
        allCircles{counter+1,5}=SVCells{i+1,4}(:,j);
        counter=counter+1;
        viscircles(SVCells{i+1,3}',SVCells{i+1,4}','Color','b');
    end
end
for i=1:size(GVCells,1)-1
    for j=1:size(GVCells{i+1,2},2)
        allCircles{counter+1,1}=counter;
        allCircles{counter+1,2}=strcat('G',num2str(GVCells{i+1,2}(1,j)));
        allCircles{counter+1,3}=strcat('G',num2str(i));
        allCircles{counter+1,4}=GVCells{i+1,3}(:,j);
        allCircles{counter+1,5}=GVCells{i+1,4}(:,j);
        counter=counter+1;
        viscircles(GVCells{i+1,3}',GVCells{i+1,4}','Color','g');
    end
end
allIntersections=getAllIntersections(allCircles,allIntersections);
[allCircles,allIntersections,allRegions]=checkAllCircles(allCircles,allIntersections,allRegions);
allRegions=checkAllIntersections(allCircles,allIntersections,allRegions);
drawRegions(allCircles,allRegions);
%disp('!!!')
end

function drawRegions(allCircles,allRegions)
eachArcStep=0.1;
for cir=2:size(allRegions,1)
    xs=[];
    ys=[];
    isNeedsDraw=false;
    for i=1:size(allRegions{cir,4},2)
        NI=allRegions{cir,4}(1,i);
        isNeedsDraw=false;
        
        if isempty(NI)==false
            xc = allCircles{NI+1,4}(1,1);
            yc = allCircles{NI+1,4}(2,1);
            r = allCircles{NI+1,5};
            
            if size(allRegions{cir,4},1)>2 && size(allRegions{cir,4},2)>2% WHY CAN THE SECOND CONDITION HAPPEN?
                if allRegions{cir,4}(5,i)==1
                    if allRegions{cir,4}(2,i)>allRegions{cir,4}(3,i)
                        theta_p = (allRegions{cir,4}(2,i))*pi/180:eachArcStep:(360)*pi/180;
                        theta_pp = (0)*pi/180:eachArcStep:(allRegions{cir,4}(3,i))*pi/180;
                        theta = [theta_p,theta_pp];
                    else
                        theta = (allRegions{cir,4}(2,i))*pi/180:eachArcStep:(allRegions{cir,4}(3,i))*pi/180;
                    end
                else
                    if allRegions{cir,4}(2,i)<allRegions{cir,4}(3,i)
                        theta_p = (allRegions{cir,4}(3,i))*pi/180:eachArcStep:(360)*pi/180;
                        theta_pp = (0)*pi/180:eachArcStep:(allRegions{cir,4}(2,i))*pi/180;
                        theta=[theta_p,theta_pp];
                    else
                        theta = (allRegions{cir,4}(3,i))*pi/180:eachArcStep:(allRegions{cir,4}(2,i))*pi/180;
                    end
                end
                
            else
                theta = (0)*pi/180:eachArcStep:(360)*pi/180;
            end
            
            if cir==5
                %disp('!!!')
            end
            
            x = xc + r*cos(theta);
            y = yc + r*sin(theta);
            if size(allRegions{cir,4},1)>2 && size(allRegions{cir,4},2)>2
                if allRegions{cir,4}(5,i)==1
                    xs=[xs,x];
                    ys=[ys,y];
                else
                    xs=[xs,flip(x)];
                    ys=[ys,flip(y)];
                    %xs=[x,xs];
                    %ys=[y,ys];
                end
            else
                xs=[xs,x];
                ys=[ys,y];
            end
            
            isNeedsDraw=true;
        end
    end
    if isNeedsDraw==true
        P = polyshape(xs,ys);
        color=[0.8,0.8,0.8];
        plot(P,'FaceColor',color,'FaceAlpha',0.6)
    end
    %disp('!!!')
end
%disp('!!!')
end

function allIntersections=getAllIntersections(allCircles,allIntersections)
counter=1;
for i=1:size(allCircles,1)-1
    for j=i+1:size(allCircles,1)-1
        [xa,yb]=circcirc(allCircles{i+1,4}(1,1),allCircles{i+1,4}(2,1),allCircles{i+1,5}(1,1),allCircles{j+1,4}(1,1),allCircles{j+1,4}(2,1),allCircles{j+1,5}(1,1));
        if ~isnan(xa(1,1))
            allIntersections{counter+1,1}=i;
            allIntersections{counter+1,2}=j;
            allIntersections{counter+1,3}=[xa(1,1);yb(1,1)];
            
            u=[allCircles{i+1,5}(1,1),0];
            v=[xa(1,1)-allCircles{i+1,4}(1,1),yb(1,1)-allCircles{i+1,4}(2,1)];
            u = u/norm(u);
            v = v/norm(v);
            angle=atan2(norm(cross([u,0],[v,0])), dot([u,0],[v,0]))*180/pi;
            if yb(1,1)-allCircles{i+1,4}(2,1)<0
                angle=360-abs(angle);
            end
            
            allIntersections{counter+1,5}=angle;
            
            u=[allCircles{j+1,5}(1,1),0];
            v=[xa(1,1)-allCircles{j+1,4}(1,1),yb(1,1)-allCircles{j+1,4}(2,1)];
            u = u/norm(u);
            v = v/norm(v);
            angle=atan2(norm(cross([u,0],[v,0])), dot([u,0],[v,0]))*180/pi;
            if yb(1,1)-allCircles{j+1,4}(2,1)<0
                angle=360-abs(angle);
            end
            
            allIntersections{counter+1,6}=angle;
            counter=counter+1;
            allIntersections{counter+1,1}=i;
            allIntersections{counter+1,2}=j;
            allIntersections{counter+1,3}=[xa(1,2);yb(1,2)];
            
            u=[allCircles{i+1,5}(1,1),0];
            v=[xa(1,2)-allCircles{i+1,4}(1,1),yb(1,2)-allCircles{i+1,4}(2,1)];
            u = u/norm(u);
            v = v/norm(v);
            angle=atan2(norm(cross([u,0],[v,0])), dot([u,0],[v,0]))*180/pi;
            if yb(1,2)-allCircles{i+1,4}(2,1)<0
                angle=360-abs(angle);
            end
            
            allIntersections{counter+1,5}=angle;
            
            u=[allCircles{j+1,5}(1,1),0];
            v=[xa(1,2)-allCircles{j+1,4}(1,1),yb(1,2)-allCircles{j+1,4}(2,1)];
            u = u/norm(u);
            v = v/norm(v);
            angle=atan2(norm(cross([u,0],[v,0])), dot([u,0],[v,0]))*180/pi;
            if yb(1,2)-allCircles{j+1,4}(2,1)<0
                angle=360-abs(angle);
            end
            
            allIntersections{counter+1,6}=angle;
            %disp('!!!')
            counter=counter+1;
        end
    end
end
%disp('!!!')
end

function allRegions=checkAllIntersections(allCircles,allIntersections,allRegions)
lastRowAllRegions=size(allRegions,1);
candidateTextNumber=1;
for i=1:size(allIntersections,1)-1
    generatedPoints=generateFourPoints(allCircles,allIntersections{i+1,3},allIntersections{i+1,1},allIntersections{i+1,2},0.05);
    %scatter(allIntersections{i+1,3}(1,1),allIntersections{i+1,3}(2,1))
    for j=1:4
        %scatter(generatedPoints(1,j),generatedPoints(2,j))
        %text(generatedPoints(1,j),generatedPoints(2,j),num2str(candidateTextNumber))
        candidateTextNumber=candidateTextNumber+1;
        forbiddens=getForbiddensOfPoint(generatedPoints(:,j),allCircles);
        forbiddenCircleIndices=[];
        if sum(forbiddens)~=0
            forbiddenCircleIndices(1,sum(forbiddens))=0;
            counter=1;
            for k=1:size(forbiddens,2)
                if forbiddens(1,k)==1
                    forbiddenCircleIndices(1,counter)=k;
                    counter=counter+1;
                end
            end
        end
        if j==1
            debug=1;
        else
            debug=0;
        end
        if sum(forbiddens)>0
            arcs=calculateRegionArc(generatedPoints(:,j),allCircles,allIntersections,forbiddens,debug);
        else
            arcs=[];
        end
        isUnique=true;
        for n=2:size(arcs,2)
            if isempty(allIntersections{arcs(4,n),4})==0
                isUnique=false;
            end
        end
        for m=1:size(allRegions,1)-1
            if isequal(allRegions{m+1,4},arcs)
                isUnique=false;
                break;
            end
        end
        
        if isUnique==true
            allRegions{lastRowAllRegions+1,3}=forbiddenCircleIndices;
            allRegions{lastRowAllRegions+1,1}=lastRowAllRegions;
            allRegions{lastRowAllRegions+1,2}=allCircles{allIntersections{i+1,1}+1,2};
            allRegions{lastRowAllRegions+1,4}=arcs;
            
            %IF THIS REGION IS UNIQUE ***
            lastRowAllRegions=lastRowAllRegions+1;
        else
            %disp('NOT UNIQUE')
        end
        
    end
    allIntersections{i+1,4}=1;
    
    %disp('!!!')
end
%disp('!!!')
end

function arcs=calculateRegionArc(internalPoint,allCircles,allIntersection,forbiddens,debug)
%scatter(internalPoint(1,1),internalPoint(2,1))
horizontalXs=[];
for i=1:size(allCircles,1)-1
    [hX,~]=linecirc(0,internalPoint(2,1),allCircles{i+1,4}(1,1),allCircles{i+1,4}(2,1),allCircles{i+1,5}(1,1));
    if ~isnan(hX(1,1))
        for m=1:size(hX,2)
            horizontalXs(1,size(horizontalXs,2)+1)=hX(1,m);
            horizontalXs(2,size(horizontalXs,2))=i;%INDEX OF THE NEIGHBOR'S CIRCLE
        end
    end
end
dists=horizontalXs(1,:)-internalPoint(1,1);
dists(dists<0)=inf;
[~,I]=min(dists);
closestInitialPoint=[horizontalXs(1,I),internalPoint(2,1)];

II=horizontalXs(2,I);

if debug==1
    %disp('!!!')
end

u=[allCircles{II+1,5}(1,1),0];
v=[closestInitialPoint(1,1)-allCircles{II+1,4}(1,1),closestInitialPoint(1,2)-allCircles{II+1,4}(2,1)];
u = u/norm(u);
v = v/norm(v);
initalAngle=atan2(norm(cross([u,0],[v,0])), dot([u,0],[v,0]))*180/pi;
if closestInitialPoint(1,2)-allCircles{II+1,4}(2,1)<0
    initalAngle=360-abs(initalAngle);
end

activeArcNI=horizontalXs(2,I);
startNeighborArc=activeArcNI;
currentActiveArc=startNeighborArc;
arcs=startNeighborArc;
currentAngle=0;

arcs(2,size(arcs,2))=initalAngle;



for i=1:size(allIntersection,1)-1
    u=[1,0];
    v=[allIntersection{i+1,3}(1,1)-internalPoint(1,1),allIntersection{i+1,3}(2,1)-internalPoint(2,1)];
    u = u/norm(u);
    v = v/norm(v);
    angle=atan2(norm(cross([u,0],[v,0])), dot([u,0],[v,0]))*180/pi;
    if allIntersection{i+1,3}(2,1)-internalPoint(2,1)<0
        angle=360-abs(angle);
    end
    allIntersection{i+1,7}=angle;
end

isNewArcFound=true;

debugCounter=0;
while(isNewArcFound==true)
    isNewArcFound=false;
    minAngle=inf;
    minAngleNotForbidden=inf;
    
    transitionArc=-1;
    transitionAngle=-1;
    transitionIntersectionIndex=-1;
    
    debugCounter=debugCounter+1;
    
    if debugCounter==4 && debug==1
        %disp('!!!!!!')
    end
    
    for g=1:size(allIntersection,1)-1
        if forbiddens(1,currentActiveArc)==1
            if allIntersection{g+1,1}(1,1)==currentActiveArc
                %if currentAngle<allIntersection{g+1,7}
                %if minAngle>allIntersection{g+1,7}
                finalAngle=-1;
                
                teta_prime=arcs(2,size(arcs,2));
                teta=allIntersection{g+1,5};
                if teta_prime>0 && teta_prime<180 && teta>0 && teta<180
                    if teta_prime<=teta
                        finalAngle=teta-teta_prime;
                    else
                        finalAngle=360-(teta_prime-teta);
                    end
                end
                if teta_prime>180 && teta_prime<360 && teta>0 && teta<180
                    finalAngle=teta+360-teta_prime;
                end
                if teta_prime>0 && teta_prime<180 && teta>180 && teta<360
                    finalAngle=-teta_prime+teta;
                end
                if teta_prime>180 && teta_prime<360 && teta>180 && teta<360
                    if teta_prime<=teta
                        finalAngle=teta-teta_prime;
                    else
                        finalAngle=360-(teta_prime-teta);
                    end
                end
                
%                 disp('From')
%                 allIntersection{g+1,1}(1,1)
%                 disp('To')
%                 allIntersection{g+1,2}(1,1)
%                 finalAngle
%                 disp('******')
                
                if debugCounter==4 && debug==1
%                     disp('!!!!!!')
                end
                if finalAngle<minAngle && finalAngle~=0% && allIntersection{g+1,7}>globalAngle
                    %candidateGlobalAngle=allIntersection{g+1,7};
                    minAngle=finalAngle;
                    transitionArc=allIntersection{g+1,2};
                    if forbiddens(1,transitionArc)==0
                        %minAngleNotForbidden=allIntersection{g+1,6};
                        %minAngle=inf;
                    else
                        %minAngleNotForbidden=inf;
                        %minAngle=allIntersection{g+1,6};
                    end
                    transitionAngle=allIntersection{g+1,6};
                    transitionIntersectionIndex=g+1;
                    arcs(3,size(arcs,2))=allIntersection{g+1,5};
                    arcs(5,size(arcs,2))=1;
                    isNewArcFound=true;
                end
                %end
                %end
            elseif allIntersection{g+1,2}==currentActiveArc
                %if currentAngle<allIntersection{g+1,7}
                %if minAngle>allIntersection{g+1,7}
                finalAngle=-1;
                
                teta_prime=arcs(2,size(arcs,2));
                teta=allIntersection{g+1,6};
                if teta_prime>0 && teta_prime<180 && teta>0 && teta<180
                    if teta_prime<=teta
                        finalAngle=teta-teta_prime;
                    else
                        finalAngle=360-(teta_prime-teta);
                    end
                end
                if teta_prime>180 && teta_prime<360 && teta>0 && teta<180
                    finalAngle=(360-teta_prime)+teta;
                end
                if teta_prime>0 && teta_prime<180 && teta>180 && teta<360
                    finalAngle=-teta_prime+teta;
                end
                if teta_prime>180 && teta_prime<360 && teta>180 && teta<360
                    if teta_prime<=teta
                        finalAngle=teta-teta_prime;
                    else
                        finalAngle=360-(teta_prime-teta);
                    end
                end
                
%                 disp('From')
%                 allIntersection{g+1,2}(1,1)
%                 disp('To')
%                 allIntersection{g+1,1}(1,1)
%                 finalAngle
%                 disp('******')
                
                if debugCounter==4 && debug==1
%                     disp('!!!!!!')
                end
                if finalAngle<minAngle && finalAngle~=0% && allIntersection{g+1,7}>globalAngle
                    %candidateGlobalAngle=allIntersection{g+1,7};
                    minAngle=finalAngle;
                    transitionArc=allIntersection{g+1,1};
                    if forbiddens(1,transitionArc)==0
                        %minAngleNotForbidden=allIntersection{g+1,5};
                        %minAngle=inf;
                    else
                        %minAngleNotForbidden=inf;
                        %minAngle=allIntersection{g+1,5};
                    end
                    transitionAngle=allIntersection{g+1,5};
                    transitionIntersectionIndex=g+1;
                    arcs(3,size(arcs,2))=allIntersection{g+1,6};
                    arcs(5,size(arcs,2))=1;
                    isNewArcFound=true;
                end
                %end
                %end
            end
        else
            if allIntersection{g+1,1}(1,1)==currentActiveArc
                %if currentAngle<allIntersection{g+1,7}
                %if minAngle>allIntersection{g+1,7}% && allIntersection{g+1,5}<minAngleNotForbidden
                finalAngle=-1;
                teta_prime=arcs(2,size(arcs,2));
                teta=allIntersection{g+1,5};
                if teta_prime>0 && teta_prime<180 && teta>0 && teta<180
                    if teta_prime>=teta
                        finalAngle=teta_prime-teta;
                    else
                        finalAngle=360-(teta-teta_prime);
                    end
                end
                if teta_prime>180 && teta_prime<360 && teta>0 && teta<180
                    finalAngle=teta_prime-teta;
                end
                if teta_prime>0 && teta_prime<180 && teta>180 && teta<360
                    finalAngle=teta_prime+360-teta;
                end
                if teta_prime>180 && teta_prime<360 && teta>180 && teta<360
                    if teta_prime>=teta
                        finalAngle=teta_prime-teta;
                    else
                        finalAngle=360-(teta-teta_prime);
                    end
                end
%                 disp('From')
%                 allIntersection{g+1,1}(1,1)
%                 disp('To')
%                 allIntersection{g+1,2}(1,1)
%                 finalAngle
%                 disp('******')
                
                if debugCounter==4 && debug==1
%                     disp('!!!!!!')
                end
                
                %angleLeftSide=min(mod(abs(minAngleNotForbidden-allIntersection{g+1,5}),360),mod(abs(360-minAngleNotForbidden+allIntersection{g+1,5}),360));
                if finalAngle<minAngleNotForbidden  && finalAngle~=0% && allIntersection{g+1,7}>globalAngle
                    %candidateGlobalAngle=allIntersection{g+1,7};
                    minAngleNotForbidden=finalAngle;
                    %minAngle=allIntersection{g+1,7};
                    transitionArc=allIntersection{g+1,2};
                    if forbiddens(1,transitionArc)==0
                        %minAngleNotForbidden=allIntersection{g+1,6};
                        %minAngle=inf;
                    else
                        %minAngleNotForbidden=inf;
                        %minAngle=allIntersection{g+1,6};
                    end
                    transitionAngle=allIntersection{g+1,6};
                    transitionIntersectionIndex=g+1;
                    arcs(3,size(arcs,2))=allIntersection{g+1,5};
                    arcs(5,size(arcs,2))=0;
                    isNewArcFound=true;
                end
                %end
                %end
            elseif allIntersection{g+1,2}==currentActiveArc
                %if currentAngle<allIntersection{g+1,7}
                %if minAngle>allIntersection{g+1,7}% && allIntersection{g+1,6}<minAngleNotForbidden
                finalAngle=-1;
                teta_prime=arcs(2,size(arcs,2));
                teta=allIntersection{g+1,6};
                if teta_prime>0 && teta_prime<180 && teta>0 && teta<180
                    if teta_prime>=teta
                        finalAngle=teta_prime-teta;
                    else
                        finalAngle=360-(teta-teta_prime);
                    end
                end
                if teta_prime>180 && teta_prime<360 && teta>0 && teta<180
                    finalAngle=teta_prime-teta;
                end
                if teta_prime>0 && teta_prime<180 && teta>180 && teta<360
                    finalAngle=teta_prime+360-teta;
                end
                if teta_prime>180 && teta_prime<360 && teta>180 && teta<360
                    if teta_prime>=teta
                        finalAngle=teta_prime-teta;
                    else
                        finalAngle=360-(teta-teta_prime);
                    end
                end
%                 disp('From')
%                 allIntersection{g+1,2}(1,1)
%                 disp('To')
%                 allIntersection{g+1,1}(1,1)
%                 finalAngle
%                 disp('******')
                
                if debugCounter==4 && debug==1
%                     disp('!!!!!!')
                end
                
                %angleLeftSide=min(mod(abs(minAngleNotForbidden-allIntersection{g+1,6}),360),mod(abs(360-minAngleNotForbidden+allIntersection{g+1,6}),360));
                if finalAngle<minAngleNotForbidden  && finalAngle~=0% && allIntersection{g+1,7}>globalAngle
                    %candidateGlobalAngle=allIntersection{g+1,7};
                    minAngleNotForbidden=finalAngle;
                    %minAngle=allIntersection{g+1,7};
                    transitionArc=allIntersection{g+1,1};
                    if forbiddens(1,transitionArc)==0
                        %minAngleNotForbidden=allIntersection{g+1,5};
                        %minAngle=inf;
                    else
                        %minAngleNotForbidden=inf;
                        %minAngle=allIntersection{g+1,5};
                    end
                    transitionAngle=allIntersection{g+1,5};
                    transitionIntersectionIndex=g+1;
                    arcs(3,size(arcs,2))=allIntersection{g+1,6};
                    arcs(5,size(arcs,2))=0;
                    isNewArcFound=true;
                end
                %end
                %end
            end
        end
    end
    if debug==1
%         disp('!!!')
    end
    if transitionArc>-1
        if arcs(3,size(arcs,2))==arcs(3,1) && size(arcs,2)>1
            break;%FOUND BEGINNING!
        end
        currentActiveArc=transitionArc;
        currentAngle=minAngle;
        %globalAngle=candidateGlobalAngle;
        arcs(1,size(arcs,2)+1)=currentActiveArc;
        arcs(2,size(arcs,2))=transitionAngle;
        arcs(4,size(arcs,2))=transitionIntersectionIndex;
        if forbiddens(1,transitionArc)==1
            minAngleNotForbidden=inf;
        else
            minAngle=inf;
        end
    end
end
arcs(3,size(arcs,2))=arcs(2,1);
% if arcs(5,size(arcs,2))==1
%     arcs(3,size(arcs,2))=arcs(2,1);
% else
%     arcs(3,size(arcs,2))=arcs(2,1);
% end
end

function outut=generateFourPoints(allCircles,intersectionPoint,circle1I,circle2I,segmentLength)
circle1Center=allCircles{circle1I+1,4};
tangent1M=((intersectionPoint(2,1)-circle1Center(2,1)))/((intersectionPoint(1,1)-circle1Center(1,1)));

% dXdY=[1;tangent1M];
% length=(dXdY(1,1)^2)+(dXdY(2,1)^2);
% dXdY=dXdY/length;
% dXdY=dXdY/100;
% point1=intersectionPoint+dXdY;
% point2=intersectionPoint-dXdY;
% line([point1(1,1),point2(1,1)],[point1(2,1),point2(2,1)])

circle1Center=allCircles{circle2I+1,4};
tangent2M=((intersectionPoint(2,1)-circle1Center(2,1)))/((intersectionPoint(1,1)-circle1Center(1,1)));

% dXdY=[1;tangent2M];
% length=(dXdY(1,1)^2)+(dXdY(2,1)^2);
% dXdY=dXdY/length;
% dXdY=dXdY/100;
% point1=intersectionPoint+dXdY;
% point2=intersectionPoint-dXdY;
% line([point1(1,1),point2(1,1)],[point1(2,1),point2(2,1)])

tangent12M=(tangent2M+tangent1M)/2;

if tangent12M>10
%     disp('!!!')
end

dXdY=[1;tangent12M];
length=(dXdY(1,1)^2)+(dXdY(2,1)^2);
dXdY=dXdY/sqrt(length);
dXdY=dXdY*segmentLength;
point1=intersectionPoint+dXdY;
point2=intersectionPoint-dXdY;
%line([point1(1,1),point2(1,1)],[point1(2,1),point2(2,1)])

tangent12M=(-1)/tangent12M;
dXdY=[1;tangent12M];
length=(dXdY(1,1)^2)+(dXdY(2,1)^2);
dXdY=dXdY/sqrt(length);
dXdY=dXdY*segmentLength;
point3=intersectionPoint+dXdY;
point4=intersectionPoint-dXdY;
%line([point3(1,1),point4(1,1)],[point3(2,1),point4(2,1)])
outut=[point1,point2,point3,point4];
end

function [allCircles,allIntersections,allRegions]=checkAllCircles(allCircles,allIntersections,allRegions)
noIntersectingCircles=1:size(allCircles,1)-1;
noIntersectingCircles(2,size(noIntersectingCircles,2))=0;
for i=1:size(allIntersections,1)-1
    noIntersectingCircles(2,allIntersections{i+1,1}(1,1))=1;
    noIntersectingCircles(2,allIntersections{i+1,2}(1,1))=1;
end
lastRowAllRegions=size(allRegions,1);
for i=1:size(noIntersectingCircles,2)
    if noIntersectingCircles(2,i)==0
        allCircles{i+1,6}=1;
        allRegions{lastRowAllRegions+1,1}=lastRowAllRegions;
        allRegions{lastRowAllRegions+1,2}=allCircles{i+1,2};
        forbiddens=getForbiddensOfPoint(allCircles{i+1,4}(:,1),allCircles);
        forbiddenCircleIndices=[];
        forbiddenCircleIndices(1,sum(forbiddens))=0;
        counter=1;
        for j=1:size(forbiddens,2)
            if forbiddens(1,j)==1
                forbiddenCircleIndices(1,counter)=j;
                counter=counter+1;
            end
        end
        allRegions{lastRowAllRegions+1,3}=forbiddenCircleIndices;
        allRegions{lastRowAllRegions+1,4}=i;
        lastRowAllRegions=size(allRegions,1);
        %disp('!!!')
    end
end
%disp('!!!')
end

function forbiddens=getForbiddensOfPoint(point,allCircles)
forbiddens(1,size(allCircles,1)-1)=0;
for i=1:size(allCircles,1)-1
    CX=allCircles{i+1,4}(1,1);
    CY=allCircles{i+1,4}(2,1);
    R=allCircles{i+1,5}(1,1);
    if (point(1,1)-CX)^2+(point(2,1)-CY)^2<R^2
        forbiddens(1,i)=1;
    end
end
end

function [SVCells,GVCells]=calcCores(shops,gyms)
dt = delaunayTriangulation(shops);
[SSV,SSR] = voronoiDiagram(dt);
dt = delaunayTriangulation(gyms);
[GSV,GSR] = voronoiDiagram(dt);

SVCells=calcNeighbors(SSV,SSR,shops,true);
GVCells=calcNeighbors(GSV,GSR,gyms,true);

SVCells=calcArcs(SVCells,shops);
GVCells=calcArcs(GVCells,gyms);
end

function VCells=calcArcs(VCells,locations)
VCells{1,5}='Arcs';
for i=1:size(VCells,1)-1
    %for j=1:size(VCells{i+1,2},2)
    %    fimplicit(@(x,y)sqrt((x-locations(VCells{i+1,2}(1,j),1))^2+(y-locations(VCells{i+1,2}(1,j),2))^2)/sqrt(((x-locations(i,1))^2+(y-locations(i,2))^2))-4)
    %end
    allIntersection=[];
    %intersectionDistances=[];
    intersections=[];
    horizontalXs=[];
    for j=1:size(VCells{i+1,2},2)
        for k=j+1:size(VCells{i+1,2},2)
            if j~=k
                %neigh1L=locations(VCells{i+1,2}(1,j),:);
                %neigh2L=locations(VCells{i+1,2}(1,k),:);
                %neigh1LX=neigh1L(1,1);
                %neigh1LY=neigh1L(1,2);
                %neigh2LX=neigh2L(1,1);
                %neigh2LY=neigh2L(1,2);
                
                %cLX=locations(i,1);
                %cLY=locations(i,2);
                
                %syms x y
                %S=solve(sqrt((x-neigh1LX)^2+(y-neigh1LY)^2)==sqrt((x-neigh2LX)^2+(y-neigh2LY)^2),y)
                %lineEqn=2*x*(neigh1LX-neigh2LX)+2*y*(neigh1LY-neigh2LY)+((neigh2LX^2)-(neigh1LX^2)+(neigh2LY^2)-(neigh1LY^2));
                %circleEqn=sqrt((x-neigh1LX)^2+(y-neigh1LY)^2)/sqrt(((x-cLX)^2+(y-cLY)^2));
                %circleEqnO=sqrt((x-neigh2LX)^2+(y-neigh2LY)^2)/sqrt(((x-cLX)^2+(y-cLY)^2));
                
                %fimplicit(y==S)
                %fimplicit(lineEqn)
                %if i==1
                %viscircles(VCells{i+1,3}',VCells{i+1,4}')
                %end
                %fimplicit(circleEqn==2)%OLD IMPLICIT WAY. NOW I HAVE CIRCLES
                %fimplicit(circleEqnO==4)
                
                %eqn=sqrt((x-locations(VCells{i+1,2}(1,j),1))^2+(y-locations(VCells{i+1,2}(1,j),2))^2)==sqrt((x-locations(VCells{i+1,2}(1,k),1))^2+(y-locations(VCells{i+1,2}(1,k),2))^2);
                %\/\/\/ OLD IMPLICIT SOLVE!
                %S=solve([lineEqn==0,circleEqn==2],[x,y]);
                %for h=1:size(S.x,1)
                %    if isreal(S.x(h,1))==true
                %        allIntersection(size(allIntersection,1)+1,1)=eval(S.x(h,1));
                %        allIntersection(size(allIntersection,1),2)=eval(S.y(h,1));
                %        intersectionDistances(size(intersectionDistances,1)+1,1)=pdist([locations(i,:);allIntersection(size(allIntersection,1),:)]);
                %    end
                %end
                %^^^ OLD IMPLICIT SOLVE!
                [xa,yb]=circcirc(VCells{i+1,3}(1,j),VCells{i+1,3}(2,j),VCells{i+1,4}(1,j),VCells{i+1,3}(1,k),VCells{i+1,3}(2,k),VCells{i+1,4}(1,k));
                if ~isnan(xa(1,1))
                    for h=1:size(xa,2)
                        allIntersection(size(allIntersection,1)+1,1)=xa(1,h);
                        allIntersection(size(allIntersection,1),2)=yb(1,h);
                        allIntersection(size(allIntersection,1),3)=VCells{i+1,2}(1,j);
                        allIntersection(size(allIntersection,1),4)=VCells{i+1,2}(1,k);
                        allIntersection(size(allIntersection,1),5)=i;
                        
                        %scatter(xa(1,h),yb(1,h));
                        %scatter(VCells{i+1,3}(1,j),VCells{i+1,3}(2,j));
                        %text(VCells{i+1,3}(1,j),VCells{i+1,3}(2,j),num2str(j))
                        
                        u=[VCells{i+1,4}(1,j),0];
                        v=[xa(1,h)-VCells{i+1,3}(1,j),yb(1,h)-VCells{i+1,3}(2,j)];
                        u = u/norm(u);
                        v = v/norm(v);
                        
                        %angle = atan2(detV,dotV)*(180/pi);
                        %aa=mod(atan2d(u,v)+360,360)
                        angle=atan2(norm(cross([u,0],[v,0])), dot([u,0],[v,0]))*180/pi;
                        if yb(1,h)-VCells{i+1,3}(2,j)<0
                            angle=360-abs(angle);
                        end
                        
                        allIntersection(size(allIntersection,1),6)=angle;
                        
                        %scatter(VCells{i+1,3}(1,k),VCells{i+1,3}(2,k));
                        %text(VCells{i+1,3}(1,k),VCells{i+1,3}(2,k),num2str(k))
                        
                        u=[VCells{i+1,4}(1,k),0];
                        v=[xa(1,h)-VCells{i+1,3}(1,k),yb(1,h)-VCells{i+1,3}(2,k)];
                        u = u/norm(u);
                        v = v/norm(v);
                        
                        %angle = atan2(detV,dotV)*(180/pi);
                        %aa=mod(atan2d(u,v)+360,360)
                        angle=atan2(norm(cross([u,0],[v,0])), dot([u,0],[v,0]))*180/pi;
                        if yb(1,h)-VCells{i+1,3}(2,k)<0
                            angle=360-abs(angle);
                        end
                        
                        allIntersection(size(allIntersection,1),7)=angle;
                        u=[VCells{i+1,4}(1,j),0];
                        v=[xa(1,h)-locations(i,1),yb(1,h)-locations(i,2)];
                        u = u/norm(u);
                        v = v/norm(v);
                        angleFromCenter=atan2(norm(cross([u,0],[v,0])), dot([u,0],[v,0]))*180/pi;
                        if yb(1,h)-locations(i,2)<0
                            angleFromCenter=360-abs(angleFromCenter);
                        end
                        allIntersection(size(allIntersection,1),8)=angleFromCenter;
                        %intersectionDistances(size(intersectionDistances,1)+1,1)=pdist([locations(i,:);allIntersection(size(allIntersection,1),[1,2])]);
                    end
                end
                %disp('!')
            end
        end
        [hX,~]=linecirc(0,locations(i,2),VCells{i+1,3}(1,j),VCells{i+1,3}(2,j),VCells{i+1,4}(1,j));
        for m=1:size(hX,2)
            horizontalXs(1,size(horizontalXs,2)+1)=hX(1,m);
            horizontalXs(2,size(horizontalXs,2))=j;%INDEX OF THE NEIGHBOR'S CIRCLE
        end
    end
    %polyout = intersect(intersections);
    %activeArcNI=0;
    dists=horizontalXs(1,:)-locations(i,1);
    dists(dists<0)=inf;
    [~,I]=min(dists);
    closestInitialPoint=[horizontalXs(1,I),locations(i,2)];
    
    II=horizontalXs(2,I);
    
    u=[VCells{i+1,4}(1,II),0];
    v=[closestInitialPoint(1,1)-VCells{i+1,3}(1,II),closestInitialPoint(1,2)-VCells{i+1,3}(2,II)];
    u = u/norm(u);
    v = v/norm(v);
    initalAngle=atan2(norm(cross([u,0],[v,0])), dot([u,0],[v,0]))*180/pi;
    if closestInitialPoint(1,2)-VCells{i+1,3}(2,II)<0
        initalAngle=360-abs(initalAngle);
    end
    
    activeArcNI=horizontalXs(2,I);
    startNeighborArc=VCells{i+1,2}(1,activeArcNI);
    currentActiveArc=startNeighborArc;
    arcs=startNeighborArc;
    currentAngle=0;
    
    arcs(2,size(arcs,2))=initalAngle;
    isNewArcFound=true;
    while(isNewArcFound==true)
        isNewArcFound=false;
        minAngle=inf;
        transitionArc=-1;
        transitionAngle=-1;
        for g=1:size(allIntersection,1)
            if allIntersection(g,3)==currentActiveArc
                if currentAngle<allIntersection(g,8)
                    if minAngle>allIntersection(g,8)
                        minAngle=allIntersection(g,8);
                        transitionArc=allIntersection(g,4);
                        transitionAngle=allIntersection(g,7);
                        arcs(3,size(arcs,2))=allIntersection(g,6);
                        isNewArcFound=true;
                    end
                end
            elseif allIntersection(g,4)==currentActiveArc
                if currentAngle<allIntersection(g,8)
                    if minAngle>allIntersection(g,8)
                        minAngle=allIntersection(g,8);
                        transitionArc=allIntersection(g,3);
                        transitionAngle=allIntersection(g,6);
                        arcs(3,size(arcs,2))=allIntersection(g,7);
                        isNewArcFound=true;
                    end
                end
            end
        end
        if transitionArc>-1
            currentActiveArc=transitionArc;
            currentAngle=minAngle;
            arcs(1,size(arcs,2)+1)=currentActiveArc;
            arcs(2,size(arcs,2))=transitionAngle;
        end
    end
    arcs(size(arcs,1),size(arcs,2))=arcs(2,1);
    VCells{i+1,5}=arcs;
    %disp('!')
end
end

function VCells=calcNeighbors(SV,SR,locations,isSelf)
VCells=cell(1,5);
VCells{1,1}='VertexIndices';
VCells{1,2}='Cell neighbors';
VCells{1,3}='Cell neighbors circle centers';
VCells{1,4}='Cell neighbors circle radiuses';
for i=1:size(SR,1)
    VCells{i+1,1}=SR{i,1};
end

for i=1:size(SR,1)
    neighbors=[];
    for j=1:size(SR{i,1},2)
        t=SR{i,1}(1,j);
        for k=1:size(SR,1)
            for m=1:size(SR{k,1},2)
                if SR{k,1}(1,m)==t
                    if k~=i && t~=1
                        neighbors(1,size(neighbors,2)+1)=k;
                    end
                end
            end
        end
    end
    VCells{i+1,2}=unique(neighbors);
end
%for i=1:size(SR,1)
%    text(locations(i,1),locations(i,2),num2str(i))
%end
%\/\/\/ DETECT NEIGHBOR CIRCLES' CENTERS AND RADIUSES
for i=1:size(SR,1)
    centers=zeros(2,size(VCells{i+1,2},2));
    radiuses=zeros(1,size(VCells{i+1,2},2));
    for j=1:size(VCells{i+1,2},2)
        centerSelf=locations(i,:);
        centerNeighbor=locations(VCells{i+1,2}(1,j),:);
        centerVector=-centerSelf+centerNeighbor;
        if isSelf==true
            centers(:,j)=(centerSelf-(1/3)*centerVector)';
        else
            centers(:,j)=(centerNeighbor+(1/3)*centerVector)';
        end
        radiuses(1,j)=sqrt((centerVector(1,1)^2)+(centerVector(1,2)^2))*(2/3);
    end
    VCells{i+1,3}=centers;
    VCells{i+1,4}=radiuses;
end
%^^^ DETECT NEIGHBOR CIRCLES' CENTERS AND RADIUSES
%disp('!')
end

function drawCore(ShopVCells,gymVCells,shops,gyms)
drawVCells(ShopVCells,'b',shops)
drawVCells(gymVCells,'g',gyms)

%disp('!!!')

end

function drawVCells(VCells,color,POIs)
h=voronoi(POIs(:,1),POIs(:,2));
set(h, 'Color', color)
eachArcStep=0.01;
for cir=2:size(VCells,1)
    xs=[];
    ys=[];
    isNeedsDraw=false;
    for i=1:size(VCells{cir,5},2)
        NI=[];
        isNeedsDraw=false;
        for j=1:size(VCells{cir,2},2)
            if VCells{cir,5}(1,i)==VCells{cir,2}(1,j)
                NI = j;
                break;
            end
        end
        if isempty(NI)==false
            xc = VCells{cir,3}(1,NI);
            yc = VCells{cir,3}(2,NI);
            r = VCells{cir,4}(1,NI);
            
            if size(VCells{cir,5},1)>2
                if VCells{cir,5}(2,i)>VCells{cir,5}(3,i)
                    theta_p = (VCells{cir,5}(2,i))*pi/180:eachArcStep:(360)*pi/180;
                    theta_pp = (0)*pi/180:eachArcStep:(VCells{cir,5}(3,i))*pi/180;
                    theta=[theta_p,theta_pp];
                else
                    theta = (VCells{cir,5}(2,i))*pi/180:eachArcStep:(VCells{cir,5}(3,i))*pi/180;
                end
            else
                theta = (0)*pi/180:eachArcStep:(360)*pi/180;
            end
            
            x = xc + r*cos(theta);
            y = yc + r*sin(theta);
            xs=[xs,x];
            ys=[ys,y];
            isNeedsDraw=true;
        end
    end
    if isNeedsDraw==true
        P = polyshape(xs,ys);
        plot(P,'FaceColor',color,'FaceAlpha',0.6)
    end
    %disp('!!!')
end
end