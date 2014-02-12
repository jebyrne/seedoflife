function [imglist] = get_imagelist(imgdir)
%--------------------------------------------------------------------------
%
% Copyright (c) 201 Jeffrey Byrne
% $Id: get_imagelist.m 79 2012-07-27 14:30:30Z jebyrne $
%
%--------------------------------------------------------------------------

%% Error checks
if (nargin ~= 1)
  error('Must include directory for list creation');
end
if (exist(imgdir, 'dir') == 0)
  error('Directory ''%s'' not found', imgdir);
end

%% For all files in the provided directory
dir_files = dir(imgdir);
num_imgs = 0;
imglist = {};
for j=1:length(dir_files)
  % Read current image filename
  [x,x,ext] = fileparts(dir_files(j).name);
  if (dir_files(j).isdir == 0) && ~isempty(ext) && ~isempty(strmatch(ext,{'.jpg','.png','.tif','.pgm','.bmp'}))
    num_imgs = num_imgs + 1;
    imglist(num_imgs) = {strcat(imgdir,filesep,dir_files(j).name)};
  end
end

% Error check
if (num_imgs == 0) 
  error(sprintf('No images found in ''%s''', imgdir));
end

