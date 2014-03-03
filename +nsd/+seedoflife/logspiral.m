function [d] = logspiral(d,di)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne
%
%--------------------------------------------------------------------------


%% Binarize - Logarithmic spiral 
D = nsd.seedoflife.reshape(d,di);
Sa = (circshift(D,[0 1 1 0]));
Sb = (circshift(D,[0 1 -1 0]));
%S = D - (0.5*Sa + 0.5*Sb);  % weighted combination
S = D - Sa;  % single log spiral
%  S(:,:,1,:) = D(:,:,1,:);  % smallest scale unnormalized (do not use circular shift back to one!)
S(:,:,1,:) = -D(:,:,1,:)+Sb(:,:,1,:);  % smallest scale reverse normalized
%  S(:,:,end,:) = D(:,:,end,:)-Sa(:,:,end,:);  % smallest scale reverse normalized
% Dm = (circshift(Dm,[0 1 1 0]));  % mask

d = reshape(S, size(d));
