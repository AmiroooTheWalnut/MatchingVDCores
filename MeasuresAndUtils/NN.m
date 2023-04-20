function d=   NN(x,y, X,Y) %calculating the distance from (x,y) to the nearest nbr at P
dmin=1000000; k=1 ;
for i=1:size(X)
    d1=(X(i)-x)^2+( Y(i)-y)^2;
    if d1 <dmin
        dmin=d1;
        k=i ;
    end
end
d =sqrt( dmin )  ;