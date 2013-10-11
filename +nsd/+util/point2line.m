function [l] = point2line(p1, p2)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
% $Id: point2line.m 91 2012-12-12 17:03:07Z jebyrne $
%
%--------------------------------------------------------------------------

%% two point form to general form
% (y-y1) = ((y2 - y1)/(x2-x1))(x-x1) 
% (x2 - x1)(y-y1) = (y2-y1)(x-x1)
% x2y - x2y1 -x1y + x1y1 = y2x - x1y2 - y1x +x1y1
% (x2-x1)y - x2y1 + x1y1 + x1y2 - x1y1 + (y1-y2)x = 0
% (x2-x1)y - x2y1 + x1y1 + x1y2 - x1y1 + (y1-y2)x = 0
% (y1-y2)x + (x2-x1)y + (x1y2 - x2y1) = 0

x1 = p1(1);
y1 = p1(2);
x2 = p2(1);
y2 = p2(2);

a = (y1-y2);
b = (x2-x1);
c = (x1*y2 - x2*y1);
 
l = [a b c]';

% l = p1 x p2 
% l = cross(p1,p2);
