function [img_xform,B] = imtransform(img,A)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne <jebyrne@cis.upenn.edu>
% $Id: imtransform.m 156 2013-04-05 13:59:15Z jebyrne $
%
%--------------------------------------------------------------------------
% xy = A*uv (xy retinal coordinates, homogeneous)
% xy = B*uv (xy pixel coordinates, homogeneous) 


%% Input
if nsd.util.iscolor(img)
  img = rgb2gray(img);  % FIXME: color interp2
end


%% Calibration
K = [1 0 size(img,2)/2; 0 1 size(img,1)/2; 0 0 1];
Kij = [1 0 size(img,1)/2; 0 1 size(img,2)/2; 0 0 1];
Binv = K*((A\eye(3))/K);
B = K*(A/K);


%% Projective transform
[M,N] = size(img);
[J,I] = meshgrid(1:N,1:M);
ij_in = [I(:) J(:)]';
xy_in = nsd.util.ij2xy(ij_in,1);
xy_out = nsd.util.dehomogenize(Binv*nsd.util.homogenize(xy_in));
ij_out = nsd.util.xy2ij(xy_out,1); 

k_valid = nsd.util.inmat(size(img),ij_out(1,:),ij_out(2,:));
Iq = ij_out(1,k_valid);
Jq = ij_out(2,k_valid);
Vq = interp2(img,Jq,Iq,'bilinear');

img_xform = zeros(size(img));
img_xform(k_valid) = Vq;


