function [imobs,imsim,imdisp,imdisc,imref,B] = observation(imrightfile,imleftfile,imdispfile,imdiscfile,imdispscale,A,opt)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne 
%
%--------------------------------------------------------------------------

%% Similarity stereo
[imobs] = nsd.preprocess(imrightfile,opt.descriptor.pp);
[imref,xx,imscale] = nsd.preprocess(imleftfile,opt.descriptor.pp);
imdisp = (1/imdispscale)*imscale*imresize(double(imread(imdispfile)), imscale, 'bicubic');  % downscaled disparity, left to right
imdisp = nsd.util.maxfilter(imdisp,3);  % boundary assignment
[imsim,B] = nsd.util.imtransform(imref,A);  % reference
imdisc = double(mat2gray(imresize(double(imread(imdiscfile)), imscale, 'bicubic')) > 0.75);  % downscaled discontinuities
  
% [A,A_prms] = nsd.util.similarity_transform([0 0], 0, 0.98); % random similarity
% imobs = nsd.util.imtransform(imobs,A);
% imsim = nsd.util.imtransform(imsim,A);

