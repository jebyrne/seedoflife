function [I] = node_index_matrix(g)
%--------------------------------------------------------------------------
%
% Copyright (c) 2009-2011 Jeffrey Byrne
%
%--------------------------------------------------------------------------

%% Input
if isstruct(g) && isfield(g,'adj')
  A = g.adj;
elseif isnumeric(g)
  A = g;
else
  error('invalid input');
end

%% Node index matrix I(k,:) = [u v] iff e_k(u,v) 
%[i,j] = find(tril(A,-1));  % lower triangular only
[i,j] = find(A);  % directed 
I = [i j];


