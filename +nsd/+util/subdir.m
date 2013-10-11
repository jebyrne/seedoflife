function [subdirlist] = subdir(indir)
%--------------------------------------------------------------------------
%
% Copyright (c) 201 Jeffrey Byrne
% $Id: subdir.m 79 2012-07-27 14:30:30Z jebyrne $
%
%--------------------------------------------------------------------------

%% Error checks
if (nargin ~= 1)
  error('Must include directory for list creation');
end
if (exist(indir, 'dir') == 0)
  error('Directory ''%s'' not found', indir);
end

%% All nontrivial subdirectories
subdir = dir(indir);
subdirlist = {};
for k=1:length(subdir)
  if subdir(k).isdir && ~any(strcmp(subdir(k).name,{'.','..'}))
    subdirlist = [subdirlist {subdir(k).name}];
  end
end

