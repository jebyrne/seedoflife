function [is] = issymmetric(A,t)
%--------------------------------------------------------------------------
%
% Copyright (c) 2009-2011 Jeffrey Byrne
% $Id: issymmetric.m 87 2012-08-28 18:42:17Z jebyrne $
%
%--------------------------------------------------------------------------

if ~exist('t','var')
  t = 1E-6;
end

is = isempty(find(abs(tril(A)' - triu(A)) > t,1));

