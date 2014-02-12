function [D] = spblkdiag(A)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne 
% $Id: spblkdiag.m 104 2013-02-04 16:45:20Z jebyrne $
%
%--------------------------------------------------------------------------


d = []; dm = []; dp = [];
for k=1:length(A)
  d = [d A{k}(1,1) A{k}(2,2)];
  dm = [dm A{k}(2,1) 0];
  dp = [dp A{k}(1,2) 0];
end 

D = spdiags([d' dm' dp'],[0 -1 1],length(d),length(d));