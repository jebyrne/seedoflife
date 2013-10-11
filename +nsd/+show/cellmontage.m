function [] = cellmontage(X)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
%
%--------------------------------------------------------------------------


%% Inputs
if ~iscell(X)
  error('input must be a cell array of images');
end


%% Cell array montage
n_img = length(X);
imgsize = size(X{1});
I = zeros(imgsize(1),imgsize(2),1,n_img);  % montage
for k=1:n_img
  I(:,:,1,k) = mat2gray(X{k});
end
montage(I);


