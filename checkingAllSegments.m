function Forbid = checkingAllSegments(gyms, shops)
% Rerurns an n by n matrix, Forbid,  where Forbid[i,j]=1 iff
%Gym i Cannot be connected to shop j. Otherwise Forbid[i,j]=0
%close all ; figure; hold on 
%plot(gyms(:,1), gyms(:,2), 'b+', 'MarkerSize', 30, 'LineWidth', 1);
%plot(shops(:,1), shops(:,2), 'r+', 'MarkerSize', 30, 'LineWidth', 1);
t(1).FontSize = 22;

% for i =1:size(gyms,1)
%   text(gyms(i,1), gyms(i,2), num2str(i))
%     text(shops(i,1), shops(i,2), num2str(i))
% 
% end

Forbid=zeros(size(gyms,1), size(shops,1) );


for i =1:size(gyms,1)
  for j=1:size(shops,1)
    Forbid(i,j)= CheckValidSegment(gyms(i,:), shops(j,:),gyms ) |...
      CheckValidSegment(shops(j,: ), gyms(i,: ) ,shops ) ;
  end
end
end

