function [imsaliency, imgrey, imblend] = demo_saliency(imfile)
%--------------------------------------------------------------------------
%
% Copyright (c) 2014 Jeffrey Byrne 
%
%--------------------------------------------------------------------------

if nargin == 0
    imfile = 'cameraman.tif';
end

imgrey = nsd.util.im2gray(imfile);
imsaliency = nsd.saliency(imfile);
imsaliency_rgb = ind2rgb(round(255*mat2gray(imsaliency)), jet(256));
imgrey_rgb = ind2rgb(round(255*imgrey), gray(256));
imblend = 0.25*imgrey_rgb + 0.75*imsaliency_rgb;

if nargout == 0
  figure(1); clf; imagesc(imblend); axis equal; axis off;
end
