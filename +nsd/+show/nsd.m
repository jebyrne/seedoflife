function [] = nsd(n_phase,do_print,outdir)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
%
%--------------------------------------------------------------------------
close all;
%do_structure = [0 1 1 1 1 1];
%do_structure = [0 0 0 0 0 1];
do_structure = [0 0 0 1 0 0 0 1];
if ~exist('n_phase','var')
  n_phase = 6;
end
if ~exist('do_print','var')
  do_print = false;
end
if ~exist('outdir','var') && do_print
  outdir = nsd.util.datefile();
end
if do_print
  nsd.util.mkdir(outdir);
  fprintf('[%s]: exporting all figures to ''%s''\n', mfilename, outdir);
end  


%% Common Parameters
n_nesting = 7;
%str_plotstyle = {'','b','g','r','c','m','y','b','g','r','c','m','y','k'};
%str_plotstyle = {'','k','r','o','y','g','b'};  
map = jet(n_nesting);  map = flipud(map);
%map = [1 0 0; 1 215/255 0; 0 100/255 0; 0 0 1; 160/255 32/255 240/255; 1 1 1];  
for k=1:n_nesting-1
  str_plotstyle{k+1} = map(k,:);
end
str_plotstyle{1} = [0.5 0.5 0.5];
str_plotstyle{n_nesting+1} = [0.5 0.5 0.5];


%% Structure 1 - Lotus 
if do_structure(1)
  figure(1); clf; hold on;  
  n_nesting = 7;
  r = 2.^[0:n_nesting-1];
  th = 0:2*pi/n_phase:(2*pi); 
  for k=length(r):-1:2
    r_phase = r(k) - r(k-1);
    [y,x] = pol2cart(th(1:end-1) + ((2*pi/n_phase)/2)*k, r_phase*ones(1,n_phase));
    for j=1:n_phase
      h=nsd.util.draw_circle([x(j) y(j)]',r(k-1),str_plotstyle{k},[],256);
      set(h,'LineWidth',1,'LineStyle','-');
      plot(x,y,'b.','Color',[0.5 0.5 0.5],'markerSize',10);
    end
  end
  for k=length(r):-1:1  
   %h=draw_circle([0 0]',r(k),'k',(r(k)/max(r))*ones(1,3));
   %h=nsd.util.draw_circle([0 0]',r(k),[0.5 0.5 0.5],[],256);
   %set(h,'LineWidth',1);
  end
  axis off; axis equal; hold off;
  if do_print
    export_fig(fullfile(outdir,sprintf('nsd_structure_01_%02d.png',n_phase)), '-transparent');
  end
end


%% Structure 2 - Bullseye
if do_structure(2)
  figure(2); clf; hold on;
  n_nesting = 7;
  r = 2.^[0:n_nesting-1];
  for k=length(r):-1:1
  %  h=draw_circle([0 0]',r(k),'k',(r(k)/max(r))*ones(1,3));
    h=nsd.util.draw_circle([0 0]',r(k),str_plotstyle{k+1},[],256);
    set(h,'LineWidth',1);
  end
  axis off; axis equal; hold off;
  if do_print
    export_fig(fullfile(outdir,sprintf('nsd_structure_02_00.png')), '-transparent');
  end
end


%% Structure 3 - Golden ratio
if do_structure(3)
  figure(3); clf; subplot(1,2,1); hold on;
  r = [0:1/64:64];
  th = 0;
  for k=2:length(r)  
    th = th + 2*pi*(((1+sqrt(5))/2)-1);
    [y(k),x(k)] = pol2cart(th,(r(k)));
    %h=nsd.util.draw_circle([x y]',1,'b',[],256);
    %set(h,'LineWidth',1);
    %round(r(k))
    %plot(x,y,'b.','Color',map(round((r(k))+1),:));
  %  pause(0.001)
  end
  plot(x,y,'b.','Color', [0 0 0]);
  axis equal; axis tight; axis off;
  drawnow;
  keyboard

  subplot(1,2,2);  hold on;
  r = 2.^[0:n_nesting/1024:n_nesting-1];
  th = 0;
  x=[]; y = [];
  for k=2:length(r)  
    th = th + 2*pi*(((1+sqrt(5))/2)-1);
    [y(k),x(k)] = pol2cart(th,(r(k)));
    h=nsd.util.draw_circle([x(k) y(k)]',log2(r(k)),'g',[],256);
    set(h,'LineWidth',1);
    %round(r(k))
    %plot(x,y,'b.','Color',map(round((r(k))+1),:));
  %  pause(0.001)
  end
  plot(x,y,'b.','Color', [0 0 0]);
  r = 2.^[0:n_nesting-1];
  for k=length(r):-1:1  
  %  h=draw_circle([0 0]',r(k),'k',(r(k)/max(r))*ones(1,3));
    h=nsd.util.draw_circle([0 0]',r(k),[0 0 1],[],256);
    set(h,'LineWidth',1);
  end
  hold off; axis equal; axis tight; axis off;
  if do_print
    export_fig(fullfile(outdir,sprintf('nsd_structure_03_00.png')), '-transparent');
  end
end




%% Structure 4 - Seed of Life
if do_structure(4)
  figure(4); clf; hold on;
  n_nesting = 7;
  %r = (((1+sqrt(5))/2)).^[0:n_nesting-1];
  %r = 1.732.^[0:n_nesting-1];
  r = 2.^[0:n_nesting-1];
  %r = [0:n_nesting-1];
  th = 0:2*pi/n_phase:(2*pi); 
  im = zeros(128,128);
  for k=length(r):-1:2
    r_phase = r(k-1);
    %[y,x] = pol2cart(th(1:end-1) + ((2*pi/n_phase)/2)*k, r_phase*ones(1,n_phase));
    [y,x] = pol2cart(th(1:end-1) + 0, r_phase*ones(1,n_phase));
    for j=1:n_phase
      h=nsd.util.draw_circle([x(j) y(j)]',r(k-1),str_plotstyle{k},[],256);
      set(h,'LineWidth',2,'LineStyle','-');
      im = im + nsd.util.imdisc(128,128,64+x(j),64+y(j),r(k-1),1,0);    
    end
  %  break
    plot(x,y,'b.','Color',[0.5 0.5 0.5],'markerSize',10);
  end
  for k=length(r):-1:1  
  %for k=[length(r) length(r)-1]  
   % h=nsd.util.draw_circle([0 0]',r(k),[0.5 0.5 0.5],[],256);
   % set(h,'LineWidth',1);
  end
  axis off; axis equal; hold off;
  if do_print
    export_fig(fullfile(outdir,sprintf('nsd_structure_04_%02d.png',n_phase)), '-transparent');
  end  
  %figure(40); imagesc(im); axis image;
end



%% Structure 5 - Clover
if do_structure(5)
  figure(5); clf; hold on;  
  n_nesting = 7;
  %r = sqrt(2).^[0:n_nesting-1];
  r = (2).^[0:n_nesting-1];
  th = 0:2*pi/n_phase:(2*pi); 
  for k=length(r)-1:-1:1
    r_phase = r(k+1);
    %[y,x] = pol2cart(th(1:end-1) + ((2*pi/n_phase)/2)*k, r_phase*ones(1,n_phase));
    [y,x] = pol2cart(th(1:end-1) + 0, r_phase*ones(1,n_phase));
    for j=1:n_phase
      h=nsd.util.draw_circle([x(j) y(j)]',r(k+1),str_plotstyle{k},[],256);
      set(h,'LineWidth',1,'LineStyle','-');
      plot(x,y,'b.','Color',[0.5 0.5 0.5],'markerSize',10);
    end
  end
  for k=length(r):-1:1  
   %h=draw_circle([0 0]',r(k),'k',(r(k)/max(r))*ones(1,3));
   %h=nsd.util.draw_circle([0 0]',r(k),[0.5 0.5 0.5],[],256);
   %set(h,'LineWidth',1);
  end
  axis off; axis equal; hold off;
  if do_print
    export_fig(fullfile(outdir,sprintf('nsd_structure_05_%02d.png',n_phase)), '-transparent');
  end
end

%% Structure 6 - Hawaiian Earring
if do_structure(6)
  figure(6); clf; hold on;
  r = (2).^[0:n_nesting-1];
  for k=1:length(r)
    h=nsd.util.draw_circle([r(k) 0]',r(k),str_plotstyle{k},[],256);
  end
  axis equal; 
  if do_print
    export_fig(fullfile(outdir,sprintf('nsd_structure_06.png')), '-transparent');
  end
end


%% Structure 7 - Logarithmic spiral
if do_structure(7)
  figure(7); clf; hold on;
  ri = (2).^[0:n_nesting-1];
  r = (2).^[0:0.01:n_nesting-1];
  a = 1;
  b = 0.8825424;  % experimentally so that diff(th) = pi/4
%  b = 0.306349;  % golden ratio
  th = log(r./a)./b;  
  [x,y] = pol2cart(th,r);
  plot(x,y,'b-');
  thi = log(ri./a)./b;  
  [x,y] = pol2cart(thi,ri);
  plot(x,y,'k.');
  axis equal; 
  diff(thi) - pi/4
  ri
  for k=1:length(ri)
    h=nsd.util.draw_circle([x(k) y(k)]',ri(k),str_plotstyle{k},[],256);
  end  
  axis off; 
  if do_print
    export_fig(fullfile(outdir,sprintf('nsd_structure_07.png')), '-transparent');
  end
end

%% Structure 8 - Logarithmic spirals
if do_structure(8)
  figure(8); clf; hold on;
  ri = (2).^[0:n_nesting-1];
  r = (2).^[0:0.01:n_nesting-1];
  a = 1;
  %b = 0.8825424;  % experimentally so that diff(th) = pi/4
  b = 0.66191;  % experimentally so that diff(th) = 2*pi/6
  th = log(r./a)./b;  
  for k=1:6
    [x,y] = pol2cart(th+((k-1)*(2*pi/6)) + (2*pi/12),r);
    plot(x,y,'g-','LineWidth',2);
    thi = ((k-1)*(2*pi/6)) + (2*pi/12) + log(ri./a)./b;
    
  
    [x,y] = pol2cart(thi,ri);
    plot(x,y,'k.');
  end
  axis equal;  axis off;
  if do_print
    export_fig(fullfile(outdir,sprintf('nsd_structure_08.png')), '-transparent');
  end
    diff(thi) - (2*pi/6)
  
end


%% NOTES
% Vesica Piscis, flower of life, lens, eye geometry, tripod of life, Metatron's Cube
% tunable to a required rotation invariance (optimal trade off) 
% use non-pooled rotations to recover precise pose


%three axes: radius, support radius, orientation

% we can use dynamic programming on SC, but the match is not what we want
% due to cost independence?  make the same argument here for optimal subproblems that we will use
% for indexing.  wrong costs, wrong distances (not perceptual distances)
% come up with counter examples.  inversions? (should be extremely
% unlikely).  out of order tower of hanoi/

%nesting is locally holistic  (data dependent weighting) 
% think accidental horse backs

% having a grid feature can be used to represent an object (non-local,
% covers whole object), but the tradeoff is that there are too many
% possible ways to find a nearest neighbor (small changes in any grid is equally important).  Use local features instead. 
% problem is that local features need larger support for improved selectivity, but as soon as we
% make them larger then we have the same object problem.  non-local. 

% think circle with a plus, circle with dash, ellipse with plus
% no hysteresis?

