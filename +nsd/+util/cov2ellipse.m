function [e] = cov2ellipse(sigma,mu,scale)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne 
% $Id: cov2ellipse.m 91 2012-12-12 17:03:07Z jebyrne $
%
%--------------------------------------------------------------------------

% e = [u v a b c]
% a(x-u)(x-u)+2b(x-u)(y-v)+c(y-v)(y-v)=1  ellipse parameterization

if ~exist('scale','var')
  scale = 1;
end

for k=1:size(sigma,2)
  Sinv = (((scale.^2)*eye(2))*reshape(sigma(:,k),[2 2]))\eye(2);
  a(k) = Sinv(1,1);
  b(k) = Sinv(1,2);
  c(k) = Sinv(2,2);
end

e = [scale*mu;a;b;c];
