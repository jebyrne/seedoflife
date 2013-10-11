function [h_graph] = graph(G,h_fig,linespec)
%--------------------------------------------------------------------------
%
% Copyright (c) 2009-2011 Jeffrey Byrne
%
%--------------------------------------------------------------------------

%% Inputs
if nargin < 2
  h_fig = figure(sum(char(mfilename))); 
end
if ~exist('linespec','var')
  linespec = 'b.-';
end
  

%% Show me!
figure(h_fig); 
hold on; axis image; 
A = G.adj + eye(size(G.adj)); % show nodes with no edges
gplot(A,G.xy,linespec);
set(gca,'YTickLabel',[],'YTick',[],'XTickLabel',[],'XTick',[]);
axis equal; axis tight; grid off;
%axis([0 max(G.xy(:,1)) 0 max(G.xy(:,2))]); axis ij;
h_graph=findobj(h_fig,'type','line');
set(h_graph,'MarkerSize',5,'LineWidth',1);
hold off;

