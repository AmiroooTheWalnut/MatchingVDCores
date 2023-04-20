function z = CheckValidSegment(A, B,  lA)
%Assumption - each segment is oriented from A_i, to B_i
%lP is a list of points, all of type A
fprintf("--\n\n\n\n\n--\n")

v = B-A ;
z=0;


for i= 1:size(lA,1)
  fprintf("i=%i ---------", i)
  C=lA(i,:);
  
  if C(1)== A(1) || C(2) ==A(2)
    fprintf("samne")
    z=0;
    continue
  end
  
  %find nearest point on the segment AB
  t= ((v * (  C-A )' )/(v*v'));
  if t<=0
    %Closest is A
    z=0;
    continue
  elseif t>1 %closest is B then
    t=1;
  end
  
  p=   A+v *  t;
  %   Detailed explanation goes here
  
  
  
  if    distance(p,A) >2 * distance(p,C)
    disp("ddddd");
    %  viscircles( C,distance(C,p) )
 %   line([p(1), C(1)], [p(2),C(2) ] ,  'LineWidth',4  );
    z=1; %violation
    return ;
    
    
  end
  
end