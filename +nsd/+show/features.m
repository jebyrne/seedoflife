function [] = features(f,h)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
%
%--------------------------------------------------------------------------

%% Inputs
if ~exist('h','var')
  h = figure(sum(char(mfilename)));
else
  figure(h);
end
  

%% Show me
set(h,'NumberTitle','off','Toolbar','none','Menubar','none','Name','Features');
% subplot(1,2,1); imagesc(f.img); colormap(gray); axis image; axis off; 
% hold on; plot(f.ij(:,2), f.ij(:,1), 'b.', 'MarkerSize', 6); hold off;
% hold on; plot(f.centroid(2), f.centroid(1), 'r+', 'MarkerSize', 20); hold off;
% title('Interest points');
imagesc(f.spyr.oe.img{1}); axis image; axis off;
title('Orientation energy');
drawnow;



