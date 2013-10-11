function [Vhat] = column_stochastic(V,e)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne <jebyrne@cis.upenn.edu>
% $Id: column_stochastic.m 203 2013-09-05 21:58:51Z jebyrne $
%
%--------------------------------------------------------------------------

%% Inputs
if nargin == 1
  e = 0;
end

%% Stochastic columns
% M = size(V,1);  % number of vectors (columns)
% N = size(V,2);  % vector dimension (rows)
% Vhat = V ./ repmat(sum(V,1)+(e.^2),M,1);  % stochastic columns
% if e == 0
%   Vhat(isnan(Vhat))=0;
% end
Vhat = scale_cols(V,1./sum(V,1));
Vhat(isnan(Vhat))=0;