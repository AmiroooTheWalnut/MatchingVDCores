clc
clear
%rng(2157)
rng(5)
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
    shopsPOIs{1,i} = images.roi.Point(gca,'Position',shops(i,:));
    addlistener(shopsPOIs{1,i},'ROIMoved',@allevents);
end
for i=1:numGyms
    gymsPOIs{1,i} = images.roi.Point(gca,'Position',gyms(i,:));
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
h=voronoi(shops(:,1),shops(:,2));
set(h, 'Color', 'r')
voronoi(gyms(:,1),gyms(:,2));
set(h, 'Color', 'g')
xlim([0,1])
ylim([0,1])
daspect([1,1,1])
numShops=size(shops,1);
for i=1:numShops
    shopsPOIs{1,i} = images.roi.Point(gca,'Position',shops(i,:));
    %addlistener(shopsPOIs{1,i},'MovingROI',@allevents);
    addlistener(shopsPOIs{1,i},'ROIMoved',@allevents);
end
numGyms=size(gyms,1);
for i=1:numGyms
    gymsPOIs{1,i} = images.roi.Point(gca,'Position',gyms(i,:));
    %addlistener(shopsPOIs{1,i},'MovingROI',@allevents);
    addlistener(gymsPOIs{1,i},'ROIMoved',@allevents);
end
[Svx,Svy]=voronoi(shops(:,1),shops(:,2));
dt = delaunayTriangulation(shops);

[SSV,SSR] = voronoiDiagram(dt);
[Svx,Svy]=voronoi(gyms(:,1),gyms(:,2));
dt = delaunayTriangulation(gyms);
[GSV,GSR] = voronoiDiagram(dt);

SVCells=calcNeighbors(SSV,SSR,shops);
GVCells=calcNeighbors(GSV,GSR,gyms);
drawCore(SVCells,GVCells);
end

function VCells=calcNeighbors(SV,SR,locations)
VCells=cell(1,5);
VCells{1,1}='VertexIndices';
VCells{1,2}='Cell neighbors';
VCells{1,3}='Cell neighbors circle centers';
VCells{1,4}='Cell neighbors circle radiuses';
VCells{1,5}='Arcs';
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
        centers(:,j)=(centerSelf-(1/3)*centerVector)';
        radiuses(1,j)=sqrt((centerVector(1,1)^2)+(centerVector(1,2)^2))*(2/3);
    end
    VCells{i+1,3}=centers;
    VCells{i+1,4}=radiuses;
end
%^^^ DETECT NEIGHBOR CIRCLES' CENTERS AND RADIUSES
for i=1:size(SR,1)
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
                neigh1L=locations(VCells{i+1,2}(1,j),:);
                neigh2L=locations(VCells{i+1,2}(1,k),:);
                neigh1LX=neigh1L(1,1);
                neigh1LY=neigh1L(1,2);
                neigh2LX=neigh2L(1,1);
                neigh2LY=neigh2L(1,2);
                
                cLX=locations(i,1);
                cLY=locations(i,2);
                
                syms x y
                %S=solve(sqrt((x-neigh1LX)^2+(y-neigh1LY)^2)==sqrt((x-neigh2LX)^2+(y-neigh2LY)^2),y)
                lineEqn=2*x*(neigh1LX-neigh2LX)+2*y*(neigh1LY-neigh2LY)+((neigh2LX^2)-(neigh1LX^2)+(neigh2LY^2)-(neigh1LY^2));
                circleEqn=sqrt((x-neigh1LX)^2+(y-neigh1LY)^2)/sqrt(((x-cLX)^2+(y-cLY)^2));
                %circleEqnO=sqrt((x-neigh2LX)^2+(y-neigh2LY)^2)/sqrt(((x-cLX)^2+(y-cLY)^2));
                
                %fimplicit(y==S)
                %fimplicit(lineEqn)
                if i==1
                    %viscircles(VCells{i+1,3}',VCells{i+1,4}')
                end
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
                        
                        dotV=u(1,1)*u(1,2)+v(1,1)*v(1,2);
                        detV=u(1,1)*v(1,2)-v(1,1)*u(1,2);
                        
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
                        
                        dotV=u(1,1)*u(1,2)+v(1,1)*v(1,2);
                        detV=u(1,1)*v(1,2)-v(1,1)*u(1,2);
                        
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
                        dotV=u(1,1)*u(1,2)+v(1,1)*v(1,2);
                        detV=u(1,1)*v(1,2)-v(1,1)*u(1,2);
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
    dotV=u(1,1)*u(1,2)+v(1,1)*v(1,2);
    detV=u(1,1)*v(1,2)-v(1,1)*u(1,2);
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
    %allIntersectionSortedFirst = sort(allIntersection,6);
    %allIntersectionSortedSecond = sort(allIntersection,7);
    %lastIntersectionIndex=0;
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
                    %                     currentAngle=allIntersectionSortedFirst(g,6);
                    %                     if allIntersectionSortedFirst(g,3)==currentActiveArc
                    %                         currentActiveArc=allIntersectionSortedFirst(g,4);
                    %                     else
                    %                         currentActiveArc=allIntersectionSortedFirst(g,3);
                    %                     end
                    %                     arcs(3,size(arcs,2))=allIntersectionSortedFirst(g,6);
                    %                     arcs(1,size(arcs,2)+1)=currentActiveArc;
                    %                     arcs(2,size(arcs,2))=allIntersectionSortedFirst(g,7);
                    %                     isNewArcFound=true;
                    %                     lastIntersectionIndex=g;
                    %                     break;
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
                    %                     currentAngle=allIntersectionSortedSecond(g,7);
                    %                     if allIntersectionSortedSecond(g,3)==currentActiveArc
                    %                         currentActiveArc=allIntersectionSortedSecond(g,4);
                    %                     else
                    %                         currentActiveArc=allIntersectionSortedSecond(g,3);
                    %                     end
                    %                     arcs(3,size(arcs,2))=allIntersectionSortedSecond(g,7);
                    %                     arcs(1,size(arcs,2)+1)=currentActiveArc;
                    %                     arcs(2,size(arcs,2))=allIntersectionSortedSecond(g,6);
                    %                     isNewArcFound=true;
                    %                     lastIntersectionIndex=g;
                    %                     break;
                end
            end
        end
        if transitionArc>-1
            currentActiveArc=transitionArc;
            currentAngle=minAngle;
            arcs(1,size(arcs,2)+1)=currentActiveArc;
            arcs(2,size(arcs,2))=transitionAngle;
        end
        %{
        for g=1:size(allIntersectionSorted,1)
            if currentAngle<allIntersectionSorted(g,activeArcOnIntersection)
                currentAngle=allIntersectionSorted(g,activeArcOnIntersection);
                if allIntersectionSorted(g,3)==currentActiveArc
                    currentActiveArc=allIntersectionSorted(g,4);
                else
                    currentActiveArc=allIntersectionSorted(g,3);
                end
                arcs(1,size(arcs,2)+1)=currentActiveArc;
                arcs(2,size(arcs,2))=allIntersectionSorted(g,activeArcOnIntersection);
                isNewArcFound=true;
                break;
            end
        end
        %}
    end
    %{
    if size(allIntersection,1)>1
        notClosest=allIntersection;
        [~,I]=min(intersectionDistances);
        notClosest(I,:)=[];
        if i==1
            scatter(notClosest(:,1),notClosest(:,2));
            scatter(allIntersection(I,1),allIntersection(I,2),'filled');
        end
    end
    %}
    arcs(size(arcs,1),size(arcs,2))=arcs(2,1);
    VCells{i+1,5}=arcs;
    %disp('!')
end
%disp('!')
end

function drawCore(ShopVCells,gymVCells)
drawVCells(ShopVCells,'r')
drawVCells(gymVCells,'g')

%disp('!!!')

end

function drawVCells(VCells,color)
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
            
            %lastAngle=(VCells{cir,5}(2,j))*pi/180;
            
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