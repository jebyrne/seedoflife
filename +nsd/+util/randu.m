function [x] = randu(rng,m,n)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
% $Id: randu.m 91 2012-12-12 17:03:07Z jebyrne $
%
%--------------------------------------------------------------------------

%% Uniform random number matrix (MxN) in the range [xmin,xmax]
x = (rng(2)-rng(1))*rand(m,n) + rng(1);
