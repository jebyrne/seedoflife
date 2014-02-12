function [A,theta,prms] = similarity_transform(tij,r,s)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne <jebyrne@cis.upenn.edu>
% $Id: similarity_transform.m 145 2013-03-11 14:10:46Z jebyrne $
%
%--------------------------------------------------------------------------


%% Default Inputs
if ~exist('s','var'), s = 1; end
if ~exist('r','var'), r  = 0; end
if ~exist('tij','var'), tij = [0 0]'; end


%% Similarity transformation matrix
R = [cos(r) -sin(r) 0; sin(r) cos(r) 0; 0 0 1];  % rotation
S = [s 0 0; 0 s 0; 0 0 1]; % scale
T = [0 0 tij(1); 0 0 tij(2); 0 0 0]; % translation
A = S*R + T; % composition


%% Similarity parameters
theta = [tij(:)',r,s];
prms.txy = tij(:);
prms.r = r;
prms.s = s;


