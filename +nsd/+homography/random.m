function [H,K] = random(M,N,s)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne 
%
%--------------------------------------------------------------------------

%% Inputs
if ~exist('s','var')
  s = 4;
end


%% Random homography 
p1 = [0 0; M 0; 0 N; M N]';
[di,dj] = pol2cart(2*pi*rand(1,4), min(M/s,N/s)*rand(1,4));
p2 = p1 + [di; dj];
p1h = nsd.util.homogenize(p1);
p2h = nsd.util.homogenize(p2);
H = nsd.homography.dlt(p1h,p2h);
K = [1 0 M/2; 0 1 N/2; 0 0 1];

%% Debugging
% figure(1); clf; hold on;
% plot(p1(2,:),p1(1,:),'r+');
% q = nsd.util.dehomogenize(H*p1h)
% plot(q(2,:),q(1,:),'g+');
% hold off;

