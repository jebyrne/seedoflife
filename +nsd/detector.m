function [fr,f] = detector(img,opt,f)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne
%
%--------------------------------------------------------------------------

%% Options
if nargin == 0
  error('invalid input');
end
if ~exist('opt','var')
  opt = nsd.opts().detector;  % defaults
end

%% Inputs
do_features = ~exist('f','var') || isempty(f);


%% Interest points
switch opt.mode    

  case {'random'}
    n = 100;
    ij = [size(img,1)*rand(1,n);size(img,2)*rand(1,n)]';
    fr = nsd.detector.frame(ij,0,1);
    
  case {'vlsift'}
    [fr_sift] = vl_sift(single(img));  
%    fr = nsd.detector.frame(round([fr_sift(2,:); fr_sift(1,:)]'),fr_sift(4,:),1);
%    fr = nsd.detector.frame(round([fr_sift(2,:); fr_sift(1,:)]'),fr_sift(4,:),1);
    fr = [fr_sift(2,:); fr_sift(1,:); fr_sift(3,:); fr_sift(4,:)];  % [i j s r]

  case {'vlsift-norot','DoG'}
    [fr_sift] = vl_sift(single(img));
    fr = [fr_sift(2,:); fr_sift(1,:); fr_sift(3,:); zeros(size(fr_sift(4,:)))];  % [i j s r]
    
  case {'HarrisLaplace','MultiscaleHarris','Hessian','MultiscaleHessian'}
     fr_covdet = vl_covdet(im2single(img),'EstimateOrientation', false,'Method',opt.mode); 
     fr = nsd.detector.frame([fr_covdet(2,:); fr_covdet(1,:)]',0,1);
     fr(4,:) = fr_covdet(3,:);
     fr(1:2,:) = round(fr(1:2,:));
     
  case 'harris'
    [xx,i,j] = nsd.util.harris(img,opt.harris.sigma,opt.harris.threshold,opt.harris.radius,0);    
    ij_interest = [i j];
    fr = nsd.detector.frame(ij_interest,0,1);        
    
  case 'canny'
    if nsd.util.iscolor(img)
      %imyuv = nsd.util.imrgb2yuv(img);
      %imedge = edge(imyuv(:,:,1),'canny') + edge(imyuv(:,:,2)+imyuv(:,:,3),'canny');  % color interest points
      [imedge,t_edge] = edge(imgrey,'canny',opt.canny.threshold);  % zerocross vs. canny?
      [i,j] = find(imedge);
      ij_interest = [i j];
    else
      %[imedge,t_edge] = edge(imgrey,'canny',opt.t_edge);  % zerocross vs. canny?
      [imedge,t_edge] = edge(img,'canny',opt.canny.threshold);  % zerocross vs. canny?
      [i,j] = find(imedge);
      ij_interest = [i j];
    end
    fr = nsd.detector.frame(ij_interest,0,1);    

  case 'canny-rot'
    if nsd.util.iscolor(img)
      %imyuv = nsd.util.imrgb2yuv(img);
      %imedge = edge(imyuv(:,:,1),'canny') + edge(imyuv(:,:,2)+imyuv(:,:,3),'canny');  % color interest points
      [imedge,t_edge] = edge(imgrey,'canny',opt.canny.threshold);  % zerocross vs. canny?
      [i,j] = find(imedge);
      ij_interest = [i j];
    else
      %[imedge,t_edge] = edge(imgrey,'canny',opt.t_edge);  % zerocross vs. canny?
      [imedge,t_edge] = edge(img,'canny',opt.canny.threshold);  % zerocross vs. canny?
      [i,j] = find(imedge);
      ij_interest = [i j];
    end
    fr = nsd.detector.frame(ij_interest,0,1);    
    fr_sift = [fr(2,:); fr(1,:); fr(3,:); fr(4,:)];  % [i j s r]
    [fr_sift] = vl_sift(im2single(img),'frames',fr_sift,'orientations');
    fr = [fr_sift(2,:); fr_sift(1,:); fr_sift(3,:); fr_sift(4,:)];  % [i j s r]
    
  case 'canny-subpixel'    
    % subpixel interest points
    error('fixme')
    imedge = bwmorph(edge(imgrey,'canny'),'dilate');  % zerocross vs. canny?
    imgreyorig = mat2gray(rgb2gray(imorig));
    imedgeorig = edge(imgreyorig,'canny');  % zerocross vs. canny?
    [i,j] = find(imedgeorig);
    scale = min(size(imgrey)) / min(size(imgreyorig));
    [k_upsmpl] = sub2ind(size(imgrey), ceil(scale*i), ceil(scale*j));  % relative to upper left, not center
    [k_upsmpl,k] = unique(k_upsmpl);
    i = i(k); j = j(k);
    k_downsmpl = find(imedge);
    [xx,xx,k] = intersect(k_downsmpl, k_upsmpl);
    ij_interest = scale*[i(k) j(k)]; % at most one interest point per pixel
    % subpixel interest points (END)
    
    % Frame
    fr = nsd.detector.frame(ij_interest,0,1);
    
  case 'harlap'
    vggdir = 'C:\jebyrne\penn\dev\nsd\deps\vgg-affine\extract_features';  % HACK!
    outfile = nsd.util.datefile();
    imwrite(img,sprintf('%s.png',outfile));
    str_exe = fullfile(vggdir,'extract_features.exe');
    str_cmd = sprintf('%s -harlap -i %s.png -sift -o1 %s_1.desc', str_exe, outfile, outfile);
    [status,result] = system(str_cmd);
    if status ~= 0 
      result
      error('invalid system command');
    end
    str_cmd = sprintf('%s -harlap -i %s.png -sift -o2 %s_2.desc', str_exe, outfile, outfile);
    [status,result] = system(str_cmd);
    if status ~= 0
      result
      error('invalid system command');
    end

    str_o1 = dlmread(sprintf('%s_1.desc', outfile));
    str_o2 = dlmread(sprintf('%s_2.desc', outfile));
    ij = str_o1(3:end,[2 1]);
    abc = str_o1(3:end,[3 4 5]);
    sr = str_o2(3:end,[4 5]);
    %fr = nsd.detector.frame(ij,0,1);    
    fr = [ij sr abc]'; 
%     imagesc(img); hold on; colormap(gray);
%     plot(fr(2,:),fr(1,:),'r+')';
        
  otherwise 
    error('undefined interest point mode ''%s''', opt.mode);
end


%% Interest points: Subsample
n_interest = size(fr,2);
if opt.n_subsample > 1
  %k_subsmpl = randperm(n_interest); % randomize
  k_subsmpl = 1:n_interest;
  k_subsmpl = k_subsmpl(1:opt.n_subsample:end); 
  fr = fr(:,k_subsmpl);
end

