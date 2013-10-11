function [fullpathname] = fullpath(pathname)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne <jebyrne@cis.upenn.edu>
% $Id: mkdir.m 79 2012-07-27 14:30:30Z jebyrne $
%
%--------------------------------------------------------------------------

% Relative path?
if pathname(1) == '.'
  fullpathname = fullfile(pwd,pathname);
else
  fullpathname = pathname;
end
