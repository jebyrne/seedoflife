function [img,imorig,s] = preprocess(filename, opt, is_forced)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne
%
%--------------------------------------------------------------------------

%% Inputs
if nargin <= 2
  is_forced = 0;  % force processing if 'filename' is numeric image
end


%% Import 
if ischar(filename)
  imorig = imread(filename);
elseif (isnumeric(filename) || islogical(filename)) && ~is_forced
  img = (filename);
  imorig = (filename);
  s = 1;
  return;  % if image input, passthrough with no modification 
elseif is_forced
  % passthrough
  imorig = filename;
  img = filename;
else
  error('invalid input');
end
if ~exist('opt','var')
  opt = nsd.opts().pp;
end
img = imorig;


%% Grayscale
if ndims(img) == 3 && opt.greyscale
  img = rgb2gray(img); 
end;


%% Image resize 
if ~isempty(opt.maxsize)
  s = opt.maxsize ./ min(size(img,1), size(img,2));  % scale factor
  img = imresize(img,s);
else
  s = 1;
end


%% Image rotate
img = imrotate(img,opt.rot);
imorig = imrotate(imorig,opt.rot);


%% Image contrast rescale [0,1]
img = mat2gray(img);
