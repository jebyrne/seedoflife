function [Vhat] = zeromean(V)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne <jebyrne@cis.upenn.edu>
% $Id: zeromean.m 79 2012-07-27 14:30:30Z jebyrne $
%
%--------------------------------------------------------------------------

%% Demean columns
M = size(V,1);  % number of vectors (columns)
N = size(V,2);  % vector dimension (rows)
Vhat = V - repmat(mean(V,1),M,1);  % zero mean columns

