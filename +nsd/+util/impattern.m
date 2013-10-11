function [im] = impattern(mode)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne 
% $Id: impattern.m 156 2013-04-05 13:59:15Z jebyrne $
%
%--------------------------------------------------------------------------


%% Patterns
switch(mode)
  case 'circle',
    im = nsd.util.imdisc(128,128,64,64,32,1,1);
  case 'square'
    im = zeros(128,128); im(32:96,32:96) = 1;
  case 'rotsquare'
    im = zeros(128,128); im(32:96,32:96) = 1;
    im = imrotate(im,16);
  case 'diamond'
    im = zeros(128,128); im(32:96,32:96) = 1;
    im = imrotate(im,45);
  case 'rectangle'
    im = zeros(128,128); im(32:96,48:80) = 1;
  case 'bar'
    im = zeros(128,128); im(:,64:end) = 1;
  case 'horizon'
    im = zeros(128,128); im(64:end,:) = 1;
  case 'line'
    im = zeros(128,128); im(64,:) = 1;
  case 'notch'
    im = zeros(128,128);
    im(64:end,:) = 1;
    im(64-16:64+16,64-32:64) = tril(ones(33,33));
    im(64-16:64+16,64:64+32) = fliplr(tril(ones(33,33)));
  case 'checkerboard'
    im = nsd.util.imsmooth(checkerboard(16));
  case 'cone'
    [I,J] = meshgrid(1:128,1:128);
    d = sqrt(sqdist([I(:) J(:)]', [64 64]'));
    im = zeros(128,128);
    im(:) = d;
  case 'nested'
    im = zeros(128,128);
    im(64-32:64+32,64-32:64+32) = 1;
    im(64-16:64+16,64-16:64+16) = 2;
    im(64-8:64+8,64-8:64+8) = 3;
    im = mat2gray(im);
  case 'blob'
    im = zeros(128,128);
    im(64,64) = 1;
    im = nsd.util.imsmooth(im,7);    
  case 'ramp'
    im = mat2gray(repmat([1:64 64*ones(1,64)],128,1));
  case 'sinusoid'
    im = mat2gray(repmat(sin(-pi:pi/64:pi),128,1));
    
  otherwise
    error('unknown pattern');
end