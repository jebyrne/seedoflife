function [A] = adjacency_radius(xy,r)
%--------------------------------------------------------------------------
%
% Copyright (c) 2009-2011 Jeffrey Byrne
%
%--------------------------------------------------------------------------

%% Neighborhood radius
n = size(xy,1);
d = sqrt(sqdist(xy',xy'));  % lightspeed
A = zeros(n,n);
for i=1:n
  % Nodes in g within radius r of i
  A(find(d(:,i) <= (r+1E-6)),i) = 1;
end
A = max(A,A');  % symmetric
A = A.*(1-eye(size(A)));  % zeros on main diagonal

