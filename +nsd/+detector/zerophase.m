function [ij] = zerophase(f, n_subpixel, t_mag, t_angle)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne
%
%--------------------------------------------------------------------------

%% Inputs
if ~exist('t_angle','var')
  t_angle = pi/8;
end
if ~exist('t_mag','var')
  t_mag = 0.1;
end
if ~exist('n_subpixel','var')
  n_subpixel = 4;
end


%% Orientation independent
b = f.spyr.b{1}{1};
for k=2:f.spyr.n_bands
   b = max(b,f.spyr.b{1}{k});
end
%b = b ./ 16;
bphase = angle(b);
bmag = abs(b);


%% Interpolation
[M,N] = size(bphase);
[yi,xi] = meshgrid(1:(1/n_subpixel):N,1:(1/n_subpixel):M);  
[y,x] = meshgrid(1:N,1:M);  
Bphase = interp2(y,x,bphase,yi,xi,'linear'); 
Bmag = interp2(y,x,bmag,yi,xi,'linear'); 


%% Thin and threshold
bwmag = double(Bmag > t_mag);
bwphase = double((Bphase.^2) < ((t_angle).^2));
skelphase = double(bwmorph(bwphase,'skel',Inf));
%skelmag = double(bwmorph(bwmag,'skel',Inf));
[k_zerocross] = find(skelphase.*bwmag);
ij = [xi(k_zerocross) yi(k_zerocross)];


%% Debugging
% figure(2); imagesc(Bphase); colorbar;
% figure(3); imagesc(Bmag); colorbar;
% figure(4); imagesc(bwphase); colormap(gray); 
% %figure(5); imagesc(imedge);
% figure(1); imagesc(f.img); colormap(gray); 
% hold on; plot(ij(:,2),ij(:,1),'r.'); hold off;
