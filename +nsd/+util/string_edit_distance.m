function [d,C] = string_edit_distance(s,t)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne 
% $Id: string_edit_distance.m 91 2012-12-12 17:03:07Z jebyrne $
%
%--------------------------------------------------------------------------
% http://www.levenshtein.net/
% string edit distance == levenshtein distance


%% Inputs
if ~ischar(s) || ~ischar(t)
  error('invalid input');
end

M = length(s)+1;
N = length(t)+1;

C = zeros(M,N);  % dynamic programming table
C(:,1) = [0:M-1];
C(1,:) = [0:N-1]';


%% Dynamic programming - bottom up
for j=2:M
 for k=2:N
    if s(j-1) == t(k-1)
      C(j,k) = C(j-1,k-1);  % match
    else
      C(j,k) = min([C(j-1,k)+1, C(j,k-1)+1, C(j-1,k-1)+1]);  % [skip s, skip t, skip both] 
    end
  end
end
d = C(end,end); % edit distance

