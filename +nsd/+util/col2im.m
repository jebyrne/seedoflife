function [img] = col2im(X,imgsize)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011-2012 Jeffrey Byrne
% $Id: col2im.m 79 2012-07-27 14:30:30Z jebyrne $
%
%--------------------------------------------------------------------------

%% Initialize
img = zeros(imgsize);
n_cols = size(X,2);
m_patch = round(sqrt(size(X,1)));


%% Zeropad
m_pad = floor(m_patch/2);
img = padarray(img,[m_pad m_pad]);
izp = m_pad+1;  % zero padding offset
jzp = m_pad+1;  % zero padding offset
ip = m_pad-1;
im = m_pad+1;
jp = m_pad-1;
jm = m_pad+1;


%% Accumulate 
[i,j] = ind2sub(imgsize,1:n_cols);
n_img = zeros(size(img));
for k=1:n_cols
  roi = reshape(X(:,k),[m_patch m_patch]);
  imin = i(k)-im+izp;
  imax = i(k)+ip+izp;
  jmin = j(k)-jm+jzp;
  jmax = j(k)+jp+jzp;
  img(imin:imax,jmin:jmax) = img(imin:imax,jmin:jmax) + roi;  % accumulate
  n_img(imin:imax,jmin:jmax) = n_img(imin:imax,jmin:jmax) + 1;
end


%% Crop
img = img ./ n_img;
img = img(m_pad+1:end-m_pad,m_pad+1:end-m_pad);  % crop padding 

