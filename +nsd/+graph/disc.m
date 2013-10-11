function [G] = disc(ij,r)
%--------------------------------------------------------------------------
%
% Copyright (c) 2009-2011 Jeffrey Byrne
%
%--------------------------------------------------------------------------

if ~exist('r','var')
  r = 1;  % unit disc graph
end

%% Outputs
xy = nsd.util.ij2xy(ij,2);
A = nsd.graph.adjacency_radius(xy,r);
G = struct('adj',A,'xy',xy,'ij',ij,'img',[]);
 
   

    
   

   

