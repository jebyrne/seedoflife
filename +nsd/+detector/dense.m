function [fr] = dense(im, dx, dy)
%--------------------------------------------------------------------------
%
% Copyright (c) 2014 Jeffrey Byrne
%
%--------------------------------------------------------------------------

[U,V] = meshgrid(dx:dx:size(im,2)-dx, dy:dy:size(im,1)-dy);
fr = [V(:) U(:)]';
fr = [fr; ones(1,size(fr,2)); zeros(1,size(fr,2))];


%% Orientation
% what if there is disagreement between scales?  multiple orientations?

%% Scale
% how to handle non-maxima scale space response?  approaching corner
% results in scale decreasing.  this is an ambiguity that needs ot be
% resolved with global matching?  regularization to bias towards unity
% scale if there is not a clear maxima? 

% degeneracies when we do not have observations across all orientations?  
% per orientation scale?  