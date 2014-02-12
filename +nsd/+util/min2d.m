function [v,ij,k,i,j] = min2d(A)
%--------------------------------------------------------------------------
%
% Copyright (c) 2009-2011 Jeffrey Byrne
% $Id: min2d.m 79 2012-07-27 14:30:30Z jebyrne $
%
%--------------------------------------------------------------------------

[v, k] = min(A(:));  % 2010a fixed this bug?
[ij] = ind2subv(size(A),k); % lightspeed toolbox

i = ij(:,1);
j = ij(:,2);

