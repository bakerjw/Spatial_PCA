function [DIST] = getDistanceMatrix(x,y)
DIST = zeros(length(x));
for i = 1:length(x)
    for j = i:length(x)
        DIST(i,j) = sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);
        DIST(j,i) = DIST(i,j);
    end
end
end