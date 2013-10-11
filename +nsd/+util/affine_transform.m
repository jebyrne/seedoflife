function [A,theta,prms] = affine_transform(txy,r,sx,sy,kx,ky)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne 
% $Id: affine_transform.m 91 2012-12-12 17:03:07Z jebyrne $
%
%--------------------------------------------------------------------------

%% Default Inputs
if ~exist('kx','var'), kx = 0; end
if ~exist('ky','var'), ky = 0; end
if ~exist('sx','var'), sx = 1; end
if ~exist('sy','var'), sy = 1; end
if ~exist('r','var'), r  = 0; end
if ~exist('txy','var'), txy = [0 0]'; end


%% Affine transformation matrix
R = [cos(r) -sin(r) 0; sin(r) cos(r) 0; 0 0 1];  % rotation
S = [sx 0 0; 0 sy 0; 0 0 1]; % scale
K = [1 ky 0; kx 1 0; 0 0 1]; % shear
T = [0 0 txy(1); 0 0 txy(2); 0 0 0]; % translation
A = K*S*R + T; % composition


%% Affine parameters
theta = [txy(:)',r,sx,sy,kx,ky];
prms.txy = txy;
prms.r = r;
prms.sx = sx;
prms.sy = sy;
prms.kx = kx;
prms.ky = ky;



