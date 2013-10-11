function [v] = indicatorvector(k,n)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
% $Id: indicatorvector.m 79 2012-07-27 14:30:30Z jebyrne $
%
%--------------------------------------------------------------------------

%% Binary indicator vector: v(1:n)=0, v(k)=1;
v = sparse(k,1,1,n,1);
