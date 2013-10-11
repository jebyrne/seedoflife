function [k_inmat] = inmat(matsize, i, j)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
% $Id: inmat.m 89 2012-09-16 23:41:37Z jebyrne $
%
%--------------------------------------------------------------------------

%% Inputs
if nargin == 2 && (size(i,2) == 2)
  ij = i;
  j = ij(:,2);
  i = ij(:,1);
end


%% Valid elements
imin = 1;
imax = matsize(1);
jmin = 1;
jmax = matsize(2);
k_inmat = find(all([(i(:) >= imin) (i(:) <= imax) (j(:) >= jmin) (j(:) <= jmax)],2));

