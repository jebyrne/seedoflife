function [h_graph] = imgraph(G,h_fig,plotstyle)
%--------------------------------------------------------------------------
%
% Copyright (c) 2009-2011 Jeffrey Byrne
%
%--------------------------------------------------------------------------

%% Inputs
if nargin < 2
  h_fig = figure(sum(char(mfilename))); 
  plotstyle = 'b-';
end

figure(h_fig); imagesc(G.img); colormap(gray); axis image; 
hold on; nsd.show.graph(G, h_fig, plotstyle); hold off;
drawnow;