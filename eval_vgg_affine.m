function [] = eval_vgg_affine(outdir, detectors, datalist)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne 
%
%--------------------------------------------------------------------------
close all;
import datasets.*;       % deps/vlbenchmarks
import benchmarks.*;     % deps/vlbenchmarks
import localFeatures.*;  % deps/vlbench



%% Outputs
if ~exist('outdir','var') || isempty(outdir)
  outdir = nsd.util.datefile();
  fprintf('[%s]: exporting to ''%s''\n', mfilename, outdir);
end
nsd.util.mkdir(outdir);


%% vlbenchmarks initialization
if ~exist('detectors','var') || isempty(detectors)
  siftDetector = localFeatures.VlFeatSift();
  detectors{1} = siftDetector;  detectors{1}.Name = 'SIFT';
end
if ~exist('datalist','var')
  datalist = {'bark', 'bikes', 'boat', 'graf', 'leuven', 'trees', 'ubc', 'wall'};  % all
end
matchBenchmark = RepeatabilityBenchmark('Mode','MatchingScore');
%matchBenchmark.Opts.overlapError = 0.8;  % JEBYRNE (TESTING) 
for k=1:length(detectors)
  detectorNames{k} = detectors{k}.Name;
end


%% vlbenchmarks run!
for j=1:length(datalist)
  dataset = VggAffineDataset('category',datalist{j});
  matchScore = zeros(numel(detectors),dataset.NumImages);
  numMatches = zeros(numel(detectors),dataset.NumImages);
  for d = 1:numel(detectors)
    for i = 2:dataset.NumImages
      [matchScore(d,i), numMatches(d,i),featCorresps,featReprojFrames] = ...
        matchBenchmark.testFeatureExtractor(detectors{d}, ...
        dataset.getTransformation(i), ...
        dataset.getImagePath(1), ...
        dataset.getImagePath(i)) ;
      
      % DEBUGGING: plot matches 
      if numMatches(d,i) > 0
        figure(4); clf;
        imshow(dataset.getImagePath(i));
        benchmarks.helpers.plotFrameMatches(featCorresps,...
          featReprojFrames,'IsReferenceImage',false,...
          'PlotMatchLine',false,'PlotUnmatched',false); drawnow;
        export_fig(fullfile(outdir,sprintf('%s_%s_%d.png',datalist{j},detectors{d}.Name,i)),'-transparent');
      end
    end
  end
  
  % Print Results
  printScores(detectorNames, matchScore*100, 'Match Score');
  printScores(detectorNames, numMatches, 'Number of matches') ;

  % Plot results
  figure(1); clf;
  plotScores(detectorNames, dataset, matchScore*100,'Matching Score'); title(''); 
  export_fig(fullfile(outdir,sprintf('%s_matchscore.png',datalist{j})),'-transparent');
  saveas(gcf,fullfile(outdir,sprintf('%s_matchscore.fig',datalist{j})));
  
  figure(2); clf;
  plotScores(detectorNames, dataset, numMatches,'Number of matches');  title(''); 
  export_fig(fullfile(outdir,sprintf('%s_nummatch.png',datalist{j})),'-transparent');
  saveas(gcf,fullfile(outdir,sprintf('%s_nummatch.fig',datalist{j})));  
  save(fullfile(outdir,sprintf('%s.mat',datalist{j})));  
end
save(fullfile(outdir,'vgg.mat'));
cd(outdir);


%% Composite score 
try, load('vgg.mat'); outdir='.'; catch ME, error('please CD into outdir to run section'); end
matchScore = zeros(length(detectors),6);
for k=1:length(datalist)
  mat = load(fullfile(outdir,sprintf('%s.mat',datalist{k})));
  matchScore = matchScore + mat.matchScore
end
compositeScore = matchScore ./ (length(datalist));
plot(1:5, 100*compositeScore(:,2:end)','+-','linewidth', 2); hold on ;
legend(detectorNames); grid on;
ylabel('Mean Matching Score');
xlabel('Perturbation');
axis([1 5 0 70]);
saveas(gcf,fullfile(outdir,'composite.fig'));
export_fig(fullfile(outdir,sprintf('composite.png')),'-transparent');
mean(compositeScore(:,2:end),2)
save(fullfile(outdir,'composite.mat'));
return;


%% Reload figures and print (optional) 
try, load('vgg.mat'); outdir='.'; catch ME, error('please CD into outdir to run section'); end
for k=1:length(datalist)
  dataset = datasets.VggAffineDataset('category',datalist{k});
  mat = load(sprintf('%s.mat',datalist{k}));
  plotScores(detectorNames, dataset, mat.matchScore*100,'Matching Score');
  %open(sprintf('%s_matchscore.fig',datalist{k}));
  title(''); 
  export_fig(sprintf('%s_matchscore.png',datalist{k}),'-transparent');
end



function plotScores(detectorNames, dataset, score, titleText)
xstart = max([find(sum(score,1) == 0, 1) + 1 1]);
xend = size(score,2);
xLabel = dataset.ImageNamesLabel;
xTicks = dataset.ImageNames;
plot(xstart:xend,score(:,xstart:xend)','+-','linewidth', 2); hold on ;
ylabel(titleText) ;
xlabel(xLabel);
set(gca,'XTick',xstart:1:xend);
set(gca,'XTickLabel',xTicks);
title(titleText);
%set(gca,'xtick',1:size(score,2));
maxScore = max([max(max(score)) 1]);
meanEndValue = mean(score(:,xend));
legendLocation = 'SouthEast';
if meanEndValue < maxScore/2
  legendLocation = 'NorthEast';
end
legend(detectorNames,'Location',legendLocation);
grid on ;
axis([xstart xend 0 maxScore]);
drawnow;


function printScores(detectorNames, scores, name)
numDetectors = numel(detectorNames);
maxNameLen = length('Method name');
for k = 1:numDetectors
  maxNameLen = max(maxNameLen,length(detectorNames{k}));
end
fprintf(['\n', name,':\n']);
formatString = ['%' sprintf('%d',maxNameLen) 's:'];
fprintf(formatString,'Method name');
for k = 2:size(scores,2)
  fprintf('\tImg#%02d',k);
end
fprintf('\n');
for k = 1:numDetectors
  fprintf(formatString,detectorNames{k});
  for l = 2:size(scores,2)
    fprintf('\t%6s',sprintf('%.2f',scores(k,l)));
  end
  fprintf('\n');
end


