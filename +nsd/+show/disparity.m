function [] = disparity(imdisp, dminmax)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
%
%--------------------------------------------------------------------------


%% Figure setup 
if ~exist('h','var')
  h = figure(sum(char(mfilename)));
end
set(h,'NumberTitle','off','Toolbar','none','Menubar','none','Name','Disparity');


%% Show me
if ~exist('dminmax','var') || isempty(dminmax)
  d = sort(nonzeros(imdisp),'descend');
  dmax = d(round(0.05*length(d)));  % max disparity is 95th percentile
  dminmax = [0 dmax];
end
imagesc(imdisp,dminmax); axis image; colorbar; axis off;
set(h,'color','white');  % for saving figure
drawnow; 

