function [imgout] = patch(img,roi,ij)
%--------------------------------------------------------------------------
%
% Copyright (c) 2009-2011 Jeffrey Byrne
% $Id: patch.m 79 2012-07-27 14:30:30Z jebyrne $
%
%--------------------------------------------------------------------------

%% Reference is center of roi
[M,N] = size(img);
[m,n] = size(roi);
imgout = padarray(img,[m n],0);  % zero padding
izp = m+1;  % zero padding offset
jzp = n+1;  % zero padding offset
ip = ceil(m/2)-1;
im = floor(m/2);
jp = ceil(n/2)-1;
jm = floor(n/2);
imin = ij(1)-im+izp;  
imax = ij(1)+ip+izp;
jmin = ij(2)-jm+jzp;  
jmax = ij(2)+jp+jzp;
imgout(imin:imax,jmin:jmax) = roi;  % insert patch
imgout = imgout(m+1:end-m,n+1:end-n);  % crop padding 

