function [A] = set_graph_adjacency(xy,r)
%--------------------------------------------------------------------------
%
% Copyright (c) 2009-2011 Jeffrey Byrne
%
%--------------------------------------------------------------------------

warning('use adjacency_radius')
A = nsd.graph.adjacency_radius(xy,r);

