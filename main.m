clc
clear
%rng(2157)
%rng(5)
%rng(7)%interesting, one circle in the middle
numShops=3;
numGyms=3;

global shops;
global gyms;

%\/\/\/ random?
shops=rand(numShops,2);
gyms=rand(numGyms,2);
%^^^ random?
%shops=[0.5,0.5;0.7,0.15;0.2,0.8;0.9,0.9];

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
figure(1)
clf
hold on
xlim([0,1])
ylim([0,1])
daspect([1,1,1])
numShops=size(shops,1);
for i=1:numShops
    shopsPOIs{1,i} = images.roi.Point(gca,'Position',shops(i,:),'Color','b');
    %addlistener(shopsPOIs{1,i},'MovingROI',@allevents);
    addlistener(shopsPOIs{1,i},'ROIMoved',@allevents);
end
numGyms=size(gyms,1);
for i=1:numGyms
    gymsPOIs{1,i} = images.roi.Point(gca,'Position',gyms(i,:),'Color','g');
    %addlistener(shopsPOIs{1,i},'MovingROI',@allevents);
    addlistener(gymsPOIs{1,i},'ROIMoved',@allevents);
end
drawForbiddenAreas(shops,gyms);
%[SVCells,GVCells]=calcCores(shops,gyms);%DONE
%drawCore(SVCells,GVCells,shops,gyms);%DONE
end

function drawForbiddenAreas(shops,gyms)
h=voronoi(shops(:,1),shops(:,2));
set(h, 'Color', 'b')

h=voronoi(gyms(:,1),gyms(:,2));
set(h, 'Color', 'g')

dt = delaunayTriangulation(shops);
[SSV,SSR] = voronoiDiagram(dt);
dt = delaunayTriangulation(gyms);
[GSV,GSR] = voronoiDiagram(dt);

SVCells=calcNeighbors(SSV,SSR,shops,false);
GVCells=calcNeighbors(GSV,GSR,gyms,false);

getAllRegions(SVCells,GVCells,SSR,GSR);

disp('!')
end

function getAllRegions(SVCells,GVCells,SSR,GSR)
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
allIntersections=cell(1,4);
allIntersections{1,1}='Circle 1 index';
allIntersections{1,2}='Circle 2 index';
allIntersections{1,3}='Point';
allIntersections{1,4}='Is checked';
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
checkAllCircles(allCircles,allIntersections,allRegions);
disp('!!!')
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
            counter=counter+1;
            allIntersections{counter+1,1}=i;
            allIntersections{counter+1,2}=j;
            allIntersections{counter+1,3}=[xa(1,2);yb(1,2)];
            counter=counter+1;
        end
    end
end
%disp('!!!')
end

function checkAllCircles(allCircles,allIntersections,allRegions)
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
        allRegions{lastRowAllRegions,1}=strcat(allCircles{i+1,2},'_F:',allCircles{i+1,3});
        allRegions{lastRowAllRegions,2}=allCircles{i+1,2};
        getForbiddensOfPoint() ...
        allRegions{lastRowAllRegions,3}=...;
        disp('!!!')
    end
end
disp('!!!')
end

function forbiddens=getForbiddensOfPoint(point,allCircles)

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
VCells=cell(1,3);
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
for i=1:size(SR,1)
    text(locations(i,1),locations(i,2),num2str(i))
end
%\/\/\/ DETECT NEIGHBOR CIRCLES' CENTERS AND RADIUSES
for i=1:size(SR,1)
    centers=zeros(2,size(VCells{i+1,2},2));
    radiuses=zeros(1,size(VCells{i+1,2},2));
    for j=1:size(VCells{i+1,2},2)
        centerSelf=locations(i,:);
        centerNeighbor=locations(VCells{i+1,2}(1,j),:);
        centerVector=-centerSelf+centerNeighbor;
        if isSelf==true
            centers(:,j)=(centerNeighbor-(1/3)*centerVector)';
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
end
end