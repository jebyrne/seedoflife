function [d] = kldivergence(mu0, sigma0, mu1, sigma1)
%--------------------------------------------------------------------------
%
% Copyright (c) 2012 Jeffrey Byrne 
% $Id: kldivergence.m 91 2012-12-12 17:03:07Z jebyrne $
%
%--------------------------------------------------------------------------


d = 0.5*(trace(sigma1\sigma0) + (mu1-mu0)'*(sigma1\eye(2))*(mu1-mu0) - log(det(sigma0)/det(sigma1)) - 2);
