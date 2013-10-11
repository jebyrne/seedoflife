function [] = run_sol_tradestudy(indir, outdir, k_study)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne 
%
%--------------------------------------------------------------------------
close all; drawnow;
do_study = zeros(1,1000);


%% Inputs
if ~exist('indir', 'var') || isempty(indir)
  indir = './data/middlebury-stereo';
end  
if ~exist('outdir','var') || isempty(outdir)
  outdir = nsd.util.datefile();
end
if ~exist('k_study','var') || isempty(k_study)
  k_study = 1:length(do_study);
end
nsd.util.mkdir(outdir);
do_study(k_study) = 1;


%% Image list
imglist = nsd.eval.simstereo.imagelist(indir);


%% Parameters
n_iter = 10;
n_deformations = 5;
y_radius = 10;
dopt = nsd.opts();
dopt.descriptor.pp.maxsize = 128;  
dopt.descriptor.features.spyr.n_scales = 5;  
dopt.correspondence.descriptor.detector.mode = 'canny';  
dopt.correspondence.descriptor.detector.n_subsample = 2;
dopt.correspondence.distance.mode = 'euclidean';
dopt.correspondence.assignment.mode = 'greedy';
dopt.correspondence.descriptor.sol.features.spyr.do_signed_orientation = false; 
dopt.correspondence.verbose = true;
dopt.correspondence.descriptor.verbose = true;


%% Test set
fprintf('[%s]: generating test set (defr=%d,iter=%d)\n', mfilename, n_deformations, n_iter);
testset = testlist(imglist, n_deformations, n_iter, pi/32, 0.5);
[n_imgs,n_deformations] = size(testset);


%% Parameters Study - orientation subbands
if do_study(10)
  fprintf('[%s]: study 1 - orientation subbands \n', mfilename);
  opt = dopt; 
  %prm = [2 4 8 16 32];
  prm = [2 3 4 6 8 16];
  eval = []; str_legend = [];
  for k=1:length(prm)
    fprintf('[%s][%d/%d]: study 1 - orientation subbands \n', mfilename, k, length(prm));
    opt.correspondence.descriptor.features.spyr.n_bands = prm(k);
    str_legend{k} = sprintf('%d orientations', prm(k));
    for j=1:n_deformations
      for i=1:n_imgs
        %[imobs,imref,imsim,imdisp,B] = testobs(testset{i,j}, opt);
        %[pr] = eval_similaritystereo(imsim, [], imobs, [], imdisp, B, y_radius, opt);
        [pr] = eval_similaritystereo(testset{i,j}, y_radius, opt);
        eval(i,j,k) = pr.y_detrate;
      end
    end
  end
  outfile = plot_similaritystereo(eval,str_legend,fullfile(outdir,'study_1_orientation_subbands.png'));
  fprintf('[%s]: study 1 - orientation subbands ''%s''\n', mfilename, outfile);
end


%% Parameter study - Lobes
if do_study(20)
  opt = dopt; 
  prm = [1 2 3 4 6 8 16];  
  eval = []; str_legend = [];
  for k=1:length(prm)
    fprintf('[%s][%d/%d]: study 2 - lobes \n', mfilename, k, length(prm));
    opt.correspondence.descriptor.sol.n_lobes = prm(k);
    str_legend{k} = sprintf('%d lobes', prm(k));
    for j=1:n_deformations
      for i=1:n_imgs
        %[imobs,imref,imsim,imdisp,B] = testobs(testset{i,j}, opt);
        %[pr] = eval_similaritystereo(imsim, [], imobs, [], imdisp, B, y_radius, opt);
        [pr] = eval_similaritystereo(testset{i,j}, y_radius, opt);
        eval(i,j,k) = pr.y_detrate;
      end
    end
  end
  outfile = plot_similaritystereo(eval,str_legend,fullfile(outdir,'study_2_lobes.png'));
  fprintf('[%s]: study 2 - lobes ''%s''\n', mfilename, outfile);  
end


%% Parameter study - Scales
if do_study(30)
  opt = dopt; 
  opt.correspondence.descriptor.detector.mode = 'canny';
  opt.correspondence.descriptor.detector.n_subsample = 16;  % TESTING  
  opt.correspondence.descriptor.pp.maxsize = 512;    
  prm = [3:7];  
  eval = []; str_legend = [];
  for k=1:length(prm)
    fprintf('[%s][%d/%d]: study 3 - scales \n', mfilename, k, length(prm));    
    opt.correspondence.descriptor.pp.maxsize = 2.^(prm(k)+2);    
    opt.correspondence.descriptor.sol.n_scales = prm(k);
    opt.correspondence.descriptor.features.spyr.n_scales = prm(k);
    str_legend{k} = sprintf('%d scales', prm(k));
    for j=1:n_deformations
      for i=1:n_imgs
        %[imobs,imref,imsim,imdisp,B] = testobs(testset{i,j}, opt);
        %[pr] = eval_similaritystereo(imsim, [], imobs, [], imdisp, B, 4*y_radius, opt);  % 4x for larger size
        [pr] = eval_similaritystereo(testset{i,j}, y_radius, opt);
        eval(i,j,k) = pr.y_detrate;
      end
    end
  end
  outfile = plot_similaritystereo(eval,str_legend,fullfile(outdir,'study_3_scales.png'));
  fprintf('[%s]: study 3 - scales ''%s''\n', mfilename, outfile);  
end



%% Parameter study - Pooling
if do_study(50)  
  opt = dopt; 
  prm = {{'max',0,'euclidean'},{'sum',0,'euclidean'},{'max',1,'euclidean'},{'sum',1,'euclidean'}};
  str_legend = {'MAX','SUM','Binary MAX','Binary SUM'};
  eval = [];
  for k=1:length(prm)
    fprintf('[%s][%d/%d]: study 5 - pooling \n', mfilename, k, length(prm));
    opt.correspondence.descriptor.features.pooling.mode = prm{k}{1};
    opt.correspondence.descriptor.sol.binarize = prm{k}{2};
    opt.correspondence.distance.mode = prm{k}{3};
    for j=1:n_deformations
      for i=1:n_imgs
        %[imobs,imref,imsim,imdisp,B] = testobs(testset{i,j}, opt);
        %[pr] = eval_similaritystereo(imsim, [], imobs, [], imdisp, B, y_radius, opt);
        [pr] = eval_similaritystereo(testset{i,j}, y_radius, opt);
        eval(i,j,k) = pr.y_detrate;
      end
    end
  end
  outfile = plot_similaritystereo(eval,str_legend,fullfile(outdir,'study_5_pooling.png'));
  fprintf('[%s]: study 5 - pooling ''%s''\n', mfilename, outfile);    
end


%% Parameter study - Signed orientations
if do_study(60)  
  opt = dopt; opt.descriptor.verbose = false; % nominal
  prm = {{1,16},{0,8}};
  str_legend = {'signed','unsigned'};
  eval = [];
  for k=1:length(prm)
    fprintf('[%s][%d/%d]: study 6 - signed orientations\n', mfilename, k, length(prm));
    opt.descriptor.features.spyr.do_signed_orientation = prm{k}{1};
    opt.descriptor.features.spyr.n_bands = prm{k}{2};
    for j=1:n_deformations
      for i=1:n_imgs
        %[imobs,imref,imsim,imdisp,B] = testobs(testset{i,j}, opt);
        %[pr] = eval_similaritystereo(imsim, [], imobs, [], imdisp, B, y_radius, opt);
        [pr] = eval_similaritystereo(testset{i,j}, y_radius, opt);
        eval(i,j,k) = pr.y_detrate;
      end
    end
  end
  outfile = plot_similaritystereo(eval,str_legend,fullfile(outdir,'study_6_signed.png'));
  fprintf('[%s]: study 6 - signed ''%s''\n', mfilename, outfile);    
end
return;



%--------------------------------------------------------------------
function [testset] = testlist(imglist, n_deformations, n_iter, maxrot, maxscale)
for k=1:length(imglist)
  % Similarity stereo
  if maxrot > 0
    rmin = [0:maxrot/n_deformations:maxrot];
  else
    rmin = zeros(1,n_deformations+1);
  end
  rmax = rmin(2:end);
  %smin = [0:0.75/n_deformations:0.75];  warning('minscale=0.25, zero rotations');
  if maxscale > 0
    smin = [0:maxscale/n_deformations:maxscale];  
  else
    smin = zeros(1,n_deformations+1);
  end
  smax = smin(2:end);
  rmin = [0 rmin];  % prepend zero deformation
  rmax = [0 rmax];  
  smin = [0 smin];
  smax = [0 smax];
  
  for j=1:(n_deformations+1)
    for i=1:n_iter
      % Geometric perturbation
      tij = [0 0];% + j*randn(1,2); 
      r = 0 + sign(randn)*((rmax(j)-rmin(j))*rand + rmin(j)); 
      s = 1 + sign(randn)*((smax(j)-smin(j))*rand + smin(j));
      
      % Reference similarity
      A = nsd.util.similarity_transform(tij,r,s);  % random similarity transform

      % Testlist
      testset{k,i,j}.imright = imglist(k).right;
      testset{k,i,j}.imleft = imglist(k).left;
      testset{k,i,j}.imdisp = imglist(k).disp;  % left to right
      testset{k,i,j}.imdisc = imglist(k).disc;
      testset{k,i,j}.imdispscale = imglist(k).dispscale;
      testset{k,i,j}.A = A;
      testset{k,i,j}.theta = [tij r s];
      testset{k,i,j}.index = [k,i,j];
    end
  end
end
testset = reshape(testset,[size(testset,1)*size(testset,2) size(testset,3)]);


%--------------------------------------------------------------------
function [pr] = eval_similaritystereo(testset,y_radius,opt,do_occlusions)

% Observation
[imright,imsimleft,imdisp,imdisc,imleft,B] = nsd.eval.simstereo.observation(testset.imright, ...
  testset.imleft, testset.imdisp, testset.imdisc, testset.imdispscale, testset.A, opt.correspondence);

% Descriptors
[ij_ref,ij_refinobs,ij_obs,w_asgn,info] = nsd.correspondence(imsimleft,imright, opt.correspondence);
fr_rgt = info.fr_obs;
fr_simlft = info.fr_ref;
k_asgn = info.k_asgn;
D = info.D;

% Ground truth
if exist('do_occlusions','var') && exist('imdisc','var') && (do_occlusions==1)
  [xx,xx,y,k_valid_simlft] = nsd.eval.simstereo.truth(fr_rgt,fr_simlft,imright,imsimleft,imdisp,imdisc,B,y_radius,1);
  keyboard
else
  [y,k_valid_simlft] = nsd.eval.simstereo.truth(fr_rgt,fr_simlft,imright,imsimleft,imdisp,[],B,y_radius,0);
end  

% Evaluation
[k_match,j_match] = intersect(k_asgn(:,1), k_valid_simlft);  % visible and matched
[k] = sub2ind(size(D),k_asgn(j_match,1),k_asgn(j_match,2));  % match linear index
D_match = D(k_match,:);  % distance matrix
y_match = y(k); % true matches
Y_match = y(k_match,:);  % truth matrix
j_correct = find(y_match == 1);
j_incorrect = find(y_match == -1);
pr.y_detrate = length(j_correct)/length(y_match);  % greedy detection rate
%[pr.recall,pr.precision,pr.info] = vl_pr(Y_match(:)',-D_match(:));  % precision-recall

% Debugging
% nsd.show.matching(imsimleft, imright, fr_simlft(1:2,k_asgn(j_match,1))', fr_rgt(1:2,k_asgn(j_match,2))'); % all matches
% nsd.show.matching(imsimleft, imright, fr_simlft(1:2,k_asgn(j_match(j_correct),1))',fr_rgt(1:2,k_asgn(j_match(j_correct),2))');  % correct matches only
% nsd.show.matching(imsimleft, imright, fr_simlft(1:2,k_asgn(j_match(j_incorrect),1))',fr_rgt(1:2,k_asgn(j_match(j_incorrect),2))');  % incorrect matches only
% drawnow

%--------------------------------------------------------------------
function [outfile] = plot_similaritystereo(eval,str_legend,outfile)
n_deformation = size(eval,2);
x = 0:(1/(n_deformation-1)):1; %x = x(2:end);% - (x(2)/2);
z_mu = squeeze(mean(eval,1));
z_var = squeeze(var(eval,0,1));
n_curves = length(str_legend);

figure(10); clf; hold on;
xlabel('Deformation'); ylabel('Detection Rate');
str_plotstyle = {'b.-','g.-','r.-','c.-','m.-','y.-','k.-','b:','g:','r:','c:','m:','y:','k:'};
for k=1:n_curves
  plot(x,z_mu(:,k),str_plotstyle{k},'LineWidth',2,'MarkerSize',15);
  %errorbar(x+0.01*k,z_mu(:,k),sqrt(z_var(:,k)),str_plotstyle{k});  % offset slightly
end
axis equal; axis([0 1 0 1]); grid on; 
hold off;
legend(str_legend);
% hline = findobj(gcf, 'type', 'line');
% set(hline,'LineWidth',2)
export_fig(outfile,'-transparent');
saveas(gcf,sprintf('%s.fig',outfile));
save(sprintf('%s.mat', outfile));


%--------------------------------------------------------------------
function [outfile] = plot_precisionrecall(recall,precision,str_legend,outfile)
str_plotstyle = {'b-','g-','r-','c-','m-','y-','k-','b:','g:','r:','c:','m:','y:','k:'};
figure(11); clf; hold on ;
n_curves = length(str_legend);
for k=1:n_curves
  plot(recall{k},precision{k},str_plotstyle{k},'linewidth',2);  
  %line([0 1], [1 1] * p / length(y), 'color', 'r', 'linestyle', '--') ;  
  %title(sprintf('precision-recall (AUC = %.2f %%)', info.auc * 100)) ;
  %legend('PR', 'random classifier', 'location', 'northwestoutside') ;
end
axis square ;
xlim([0 1]) ; xlabel('recall') ;
ylim([0 1]) ; ylabel('precision') ;
legend(str_legend,'location', 'northwestoutside');
grid on; drawnow;
export_fig(outfile,'-transparent');
hold off;
saveas(gcf,sprintf('%s.fig',outfile));
save(sprintf('%s.mat', outfile));
return;



%% Print me
figname = dir('*.fig'); figname = {figname.name};
for k=1:length(figname)
  open(figname{k});
  hline = findobj(gcf, 'type', 'line');
  set(hline,'LineWidth',2)
  export_fig(sprintf('%s.png',figname{k}),'-transparent');
end


