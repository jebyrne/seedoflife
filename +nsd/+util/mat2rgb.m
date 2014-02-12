function [rgb] = mat2rgb(A,cmap)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
% $Id: mat2rgb.m 91 2012-12-12 17:03:07Z jebyrne $
%
%--------------------------------------------------------------------------

if size(A,3) ~= size(cmap,1)
  error('invalid input');
end
if nargin < 2
  cmap = jet(size(A,3));
end

[a_max,j] = max(A,[],3);  % display only
rgb = ind2rgb(reshape(j,[size(A,1) size(A,2)]),cmap);
rgb = repmat(mat2gray(a_max),[1 1 3]).*rgb;
