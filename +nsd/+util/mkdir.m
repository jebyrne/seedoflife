function [] = mkdir(dirname)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne <jebyrne@cis.upenn.edu>
% $Id: mkdir.m 79 2012-07-27 14:30:30Z jebyrne $
%
%--------------------------------------------------------------------------

if ~exist(dirname,'dir')
  mkdir(dirname);
end