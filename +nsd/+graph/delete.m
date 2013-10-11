function [g] = delete(g,k_node,k_edge)
%--------------------------------------------------------------------------
%
% Copyright (c) 2009-2011 Jeffrey Byrne
%
%--------------------------------------------------------------------------

k = setdiff(1:nsd.graph.size(g), k_node);  % nodes remaining
g.adj = g.adj(k,k);
g.xy = g.xy(k,:);
g.ij = g.ij(k,:);
g.W = g.W(k,k);