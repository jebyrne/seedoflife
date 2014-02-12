function [v,ij,k,i,j] = max2d(A)
%--------------------------------------------------------------------------
%
% Copyright (c) 2009-2011 Jeffrey Byrne
% $Id: max2d.m 79 2012-07-27 14:30:30Z jebyrne $
%
%--------------------------------------------------------------------------

%% Faster
[v, k] = max(A(:));  
[ij] = ind2subv(size(A),k); % lightspeed toolbox
i = ij(:,1);
j = ij(:,2);



%% Slower
% %% Two dimensional max
% [a,ai] = max(A,[],1);
% [v,aj] = max(a);
% 
% %% Outputs
% i=ai(aj);
% j=aj;


