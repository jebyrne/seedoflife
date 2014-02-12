function [Vhat,Nhat] = normalize(V,e)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne <jebyrne@cis.upenn.edu>
% $Id: normalize.m 117 2013-02-13 02:15:12Z jebyrne $
%
%--------------------------------------------------------------------------

%% Inputs
if nargin == 1
  e = 0;
end

%% Robustified norm
% v = v / sqrt(|v|^2 + e^2)
M = size(V,1);  % vector dimension 
N = size(V,2);  % number of vectors 
%Vhat = V ./ repmat(sqrt(sum(V.^2,1)+(e.^2)),M,1);  % unit norm columns
Nhat = sqrt(sum(V.^2,1)+(e.^2));  % normalization factor
Vhat = scale_cols(V, 1./Nhat);  % unit norm columns (faster)
if e == 0
  Vhat(isnan(Vhat))=0;
end

