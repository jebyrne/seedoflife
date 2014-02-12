function [G] = new(ij, A, img)
%--------------------------------------------------------------------------
%
% Copyright (c) 2009-2011 Jeffrey Byrne
%
%--------------------------------------------------------------------------

n = size(ij,1);
xy = fliplr(ij);
if nargin == 1
  A = sparse(n,n);
  img = [];
end
G = struct('adj',A,'xy',xy,'ij',ij,'W',A,'img',img);
 
   

    
   

   

