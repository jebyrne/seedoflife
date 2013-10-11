function [k,k_valid] = sub2ind(matsize,i,j)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011-2012 Jeffrey Byrne
% $Id: sub2ind.m 156 2013-04-05 13:59:15Z jebyrne $
%
%--------------------------------------------------------------------------

%% Inputs
if nargin == 2 && (size(i,2) == 2)
  ij = i;
  j = ij(:,2);
  i = ij(:,1);
end


%% Valid?
%k_valid = nsd.util.inrect([i j], [1 1 matsize(2) matsize(1)]);
k_valid = nsd.util.inmat(matsize,i,j);


%% Index
if ~isempty(k_valid)
  n_index = length(i);
  k = NaN(n_index,1);
  %k(k_valid) = sub2ind(matsize,i(k_valid),j(k_valid));
  k = subv2ind(matsize(1:2),[nsd.util.columnize(i(k_valid)),nsd.util.columnize(j(k_valid))]);
else
  k = [];
end
