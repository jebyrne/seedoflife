function [G] = initialize(ij,r,img)
%--------------------------------------------------------------------------
%
% Copyright (c) 2009-2011 Jeffrey Byrne
%
%--------------------------------------------------------------------------

%% Inputs
if ~exist('img','var') || isempty(img);
  img = sparse(ceil(max(ij(:,1))),ceil(max(ij(:,2))));
end

%% Outputs
xy = fliplr(ij);
A = nsd.graph.adjacency_radius(xy,r);
G = struct('adj',A,'xy',xy,'ij',ij,'img',img,'W',A);
[G.n_nodes,G.n_edges] = nsd.graph.size(G); 
   

    
   

   

