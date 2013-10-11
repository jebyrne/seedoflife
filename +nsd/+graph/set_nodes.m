function [G] = set_nodes(ij)
%--------------------------------------------------------------------------
%
% Copyright (c) 2009-2011 Jeffrey Byrne
%
%--------------------------------------------------------------------------

n = size(ij,1);
xy = fliplr(ij);
A = sparse(n,n);
G = struct('adj',A,'xy',xy,'ij',f.ij);
 
   

    
   

   

