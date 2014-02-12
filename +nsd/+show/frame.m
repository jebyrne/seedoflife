function [] = frame(img,fr,h)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
%
%--------------------------------------------------------------------------

%% Inputs
if ~exist('linewidth','var')
  linewidth = 1;
end


%% Figure setup
if ~exist('h','var')
  h = figure(sum(char(mfilename)));
  set(h,'NumberTitle','off','Toolbar','none','Menubar','none','Name','Frame');
end


%% Show me
imagesc(img); colormap(gray); axis image; axis off;
hold on;
plot(fr(2,:),fr(1,:),'g.');
hold off;


%% Show me (vlfeat)
fr_sift = [fr(2,:); fr(1,:); fr(3,:); fr(4,:)];
hold on;  h1 = vl_plotframe(fr_sift); hold off;
set(h1,'color','g','linewidth',3);

