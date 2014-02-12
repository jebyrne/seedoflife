function [g] = subgraph(g,k_subgraph)
%--------------------------------------------------------------------------
%
% Copyright (c) 2009-2011 Jeffrey Byrne
%
%--------------------------------------------------------------------------

k_nbrs = nsd.graph.neighbors(g,k_subgraph);
k_subgraph = union(k_nbrs, k_subgraph);
g = nsd.graph.delete(g, setdiff(1:nsd.graph.size(g), k_subgraph));
