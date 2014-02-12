function [H] = dlt(p1h, p2h)
%--------------------------------------------------------------------
%
% Author: Jeffrey Byrne (jbyrne@ssci.com)
%
%--------------------------------------------------------------------


%% Input check
if (size(p1h, 1) == 2)
  p1h = nsd.util.homogenize(p1h);
end
if (size(p2h, 1) == 2)
  p2h = nsd.util.homogenize(p2h);
end
if size(p1h, 2) ~= size(p2h, 2)
  error('[dlt]: Inconsistent point correspondence');
end


%% Normalization
mu = mean(p1h(1:2,:),2);
K1 = [1 0 0; 0 1 0; -mu' 1]';  % translation
s = mean(sqrt(sum((nsd.util.dehomogenize(K1*p1h)).^2,1)));  % mean deviation 
K1 = [sqrt(2)/s 0 0; 0 sqrt(2)/s 0; -(sqrt(2)/s)*mu' 1]';  % conditioning
p1hn = K1*p1h;  % normalization

mu = mean(p2h(1:2,:),2);
K2 = [1 0 0; 0 1 0; -mu' 1]';  % translation
s = mean(sqrt(sum((nsd.util.dehomogenize(K2*p2h)).^2,1)));  % mean deviation 
K2 = [sqrt(2)/s 0 0; 0 sqrt(2)/s 0; -(sqrt(2)/s)*mu' 1]';  % conditioning
p2hn = K2*p2h;  % normalization


%% Direct Linear Transform
A = [];
for k=1:size(p1hn,2)
  X = p1hn(:,k)'; x=p2hn(1,k); y=p2hn(2,k); w=p2hn(3,k);
  Ai = [0 0 0 -w.*X y.*X; w.*X 0 0 0 -x.*X; -y.*X x.*X 0 0 0];
  A = [A; Ai];
end


%% TESTING COVARIANCE
% for k=1:(size(A,1)/3)
%   sigma{k} = eye(3);
%   sigma{k}(1:2,1:2) = rand(2,2) + eye(2);
% end
% S = chol(blkdiag(sigma{:}),'lower');
% A = S\A;

% DLT!
[U,S,V] = svd(A);
H = reshape(V(:,9), 3, 3)';  % unknown matrix elements are vectorized rowwise

% Choose the correct sign of H such that y'Hx > 0
if (p2hn(:,1)'*H*p1hn(:,1) < 0)
  H = -H;
end

% Correction 
H = K2\(H*K1);

% Normalize by second largest singular value
[U,S,V] = svd(H);
H = H ./ S(2,2);


