function [x] = line_line_intersection(p1, q1, p2, q2)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne 
% $Id: line_line_intersection.m 159 2013-04-17 21:00:53Z jebyrne $
%
%-------------------------------------------------------------------------- 
% Given lines l1 defined by two points (p1,q1) and line l2 defined by two 
% points (p2,q2) return the point of intersection or NaN if the two lines 
% do not intersect


%% Line-line intersection
% point coordinates (1)
x1 = p1(1);
y1 = p1(2);
x2 = q1(1);
y2 = q1(2);

% point coordinates (n)
x3 = p2(1,:);
y3 = p2(2,:);
x4 = q2(1,:);
y4 = q2(2,:);

% Analytic solution (method of determinants)
x(1,:) = ((x1.*y2 - y1.*x2)*(x3-x4) - (x1-x2).*(x3.*y4-y3.*x4))./((x1-x2).*(y3-y4)-(y1-y2).*(x3-x4));
x(2,:) = ((x1.*y2 - y1.*x2)*(y3-y4) - (y1-y2).*(x3.*y4-y3.*x4))./((x1-x2).*(y3-y4)-(y1-y2).*(x3-x4));

%if any(~isfinite(x))
% if (x(1)==NaN) || (x(2)==NaN)
%   x = NaN(2,1);
% end

