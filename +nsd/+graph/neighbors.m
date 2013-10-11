function [k_nbrs, n_nbrs] = neighbors(g,k)
%--------------------------------------------------------------------------
%
% Copyright (c) 2009-2011 Jeffrey Byrne
%
%--------------------------------------------------------------------------

[xx,k_nbrs] = find(g.adj(k,:)); 
k_nbrs = unique(k_nbrs(:));
n_nbrs = length(k_nbrs);   

