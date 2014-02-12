function [fr] = frame(ij, r, s)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne
%
%--------------------------------------------------------------------------

%% Inputs
n_interestpoint = size(ij,1);
if isscalar(r)
  r = r*ones(n_interestpoint,1);
end
if isscalar(s)
  s = s*ones(n_interestpoint,1);
end
  

%% Frame
fr = [ij s(:) r(:)]';

