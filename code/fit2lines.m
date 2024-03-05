function [a, r, sel] = fit2lines(x, y)

% temporarily rotate the data, to create a 'cross'
[u, s, v] = svd([x y], 'econ');

x_ = [x y]*v(:,1);
y_ = [x y]*v(:,2);

% first guess is fit a line to the points in the left-lower/right-upper
% quadrant, and the left-upper/right-lower quadrant

sel = (x_<0 & y_<0) | (x_>0 & y_>0);

p1  = polyfit(x(sel),  y(sel),  1);
p2  = polyfit(x(~sel), y(~sel), 1);

r = (sum((y(sel)-x(sel).*p1(1)-p1(2)).^2)+sum((y(~sel)-x(~sel).*p2(1)-p2(2)).^2))./(y'*y);

a = [p1(1) p2(1)];
