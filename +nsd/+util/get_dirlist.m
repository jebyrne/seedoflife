function [dirlist] = get_dirlist(indir)
%--------------------------------------------------------------------------
%
% Copyright (c) 201 Jeffrey Byrne
% $Id: get_dirlist.m 79 2012-07-27 14:30:30Z jebyrne $
%
%--------------------------------------------------------------------------

%% Error checks
if (nargin ~= 1)
  error('Must include directory for list creation');
end
if (exist(indir, 'dir') == 0)
  error('Directory ''%s'' not found', indir);
end

%% For all files in the provided directory
dir_files = dir(indir);
num_dirs = 0;
dirlist = {};
for j=1:length(dir_files)
  if (dir_files(j).isdir == 1) && ~strcmp(dir_files(j).name,'.') && ~strcmp(dir_files(j).name,'..')
    num_dirs = num_dirs + 1;
    dirlist(num_dirs) = {dir_files(j).name};
  end
end


