function [fr] = dense(im, dx, dy)
%--------------------------------------------------------------------------
%
% Copyright (c) 2014 Jeffrey Byrne
%
%--------------------------------------------------------------------------

[U,V] = meshgrid(1:dx:size(im,2), 1:dy:size(im,1));
fr = [V(:) U(:)]';
fr = [fr; ones(1,size(fr,2)); zeros(1,size(fr,2))];
