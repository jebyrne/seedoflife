function [k_closest] = nearest_interestpoint(im,fr,h)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
%
%--------------------------------------------------------------------------

%% Figure setup
if ~exist('h','var')
  h = figure(sum(char(mfilename)));
  imagesc(im); axis image;  axis off;  colormap(gray);
  hold on; plot(fr(2,:),fr(1,:),'b.'); hold off;  
  h=gca;
end
axes(h);


%% Closest interestpoint to selection
fprintf('[nsd.show.%s]: click on an interest point in the image\n', mfilename);
[j,i] = ginput(1);
ij_sel = [i j];
d = sqdist(fr(1:2,:), ij_sel');
[y,k_closest] = min(d);


%% Plot closest
hold on; plot(fr(2,k_closest),fr(1,k_closest),'r+'); hold off;
title(sprintf('interest point = %d', k_closest));






