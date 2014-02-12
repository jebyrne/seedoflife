function [imgout] = imcrop(img,rect)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
% $Id: imcrop.m 91 2012-12-12 17:03:07Z jebyrne $
%
%--------------------------------------------------------------------------
% imcrop with symmetric padding for rect outside img bounds
% rect = [x y width height]

[N,M] = size(img);
if min(rect(1:2)) >= 1 && (rect(3) <= M) && (rect(4) <= N)
  % within image - just use image processing toolbox
  imgout = imcrop(img,rect);
else 
  % zero pad crop 
  imgout = imcrop(img,rect); 
  if rect(1) < 1    
    imgout = padarray(imgout,[0 -rect(1)],'replicate','pre');
  end
  if rect(2) < 1
    imgout = padarray(imgout,[-rect(2) 0],'replicate','pre');
  end
  if rect(1)+rect(3) > M
    imgout = padarray(imgout,[0 (rect(1)+rect(3)-M)],'replicate','post');
  end
  if rect(2)+rect(4) > N
    imgout = padarray(imgout,[(rect(2)+rect(4)-N) 0],'replicate','post');
  end
end
