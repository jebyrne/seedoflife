function [V,E] = size(g)
%--------------------------------------------------------------------------
%
% Copyright (c) 2009-2011 Jeffrey Byrne
%
%--------------------------------------------------------------------------

V = size(g.adj,1);  % number of nodes
%E = nnz(tril(g.adj,-1)); % number of edges (not including self edges)
E = nnz(g.adj); % number of directed edges 

