function [is] = isfile(filename)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne 
% $Id: isfile.m 79 2012-07-27 14:30:30Z jebyrne $
%
%--------------------------------------------------------------------------

if ischar(filename) && ~(exist(filename,'file') == 0) 
  is = true;
else
  is = false;
end