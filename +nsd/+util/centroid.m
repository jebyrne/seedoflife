function [ij,k] = centroid(A)
%--------------------------------------------------------------------------
%
% Copyright (c) 2009-2011 Jeffrey Byrne
% $Id: centroid.m 79 2012-07-27 14:30:30Z jebyrne $
%
%--------------------------------------------------------------------------

%% Subscript
[j,i] = meshgrid(1:size(A,2),1:size(A,1));
a = A(:) ./ sum(A(:));
ij(1) = i(:)'*a;
ij(2) = j(:)'*a;

%% Index
ij_round = round(ij);
ij_round(1) = max(min(ij_round(1),size(A,1)),1);  % crop
ij_round(2) = max(min(ij_round(2),size(A,2)),1);  % crop
try
  k = sub2ind(size(A),ij_round(1), ij_round(2));
catch ME
  keyboard
end
