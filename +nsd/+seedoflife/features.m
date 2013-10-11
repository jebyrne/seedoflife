function [f, imgrey] = features(img,opt)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne
%
%--------------------------------------------------------------------------

%% Options
if ~exist('opt','var')
  opt = nsd.opts().features;  % defaults
end
if ~exist('img','var')
  imorig = imread('peppers.png');  % test image
  img = nsd.preprocess(imorig,opt.dmin);  % test image
end
if nsd.util.isfile(img)
  imorig = imread(img);
  %imorig = sub_imread(img);  warning('using rvc_imread');
  img = nsd.preprocess(imorig, opt.dmin, 0);
elseif ischar(img)
  error('image ''%s'' not found', img)
else
  imorig = img;
end


%% Greyscale
if nsd.util.iscolor(img)
  imgrey = mat2gray(rgb2gray(img));
else
  imgrey = img;
end


%% Steerable pyramid
spyr = nsd.seedoflife.cspyr(img, opt.spyr);  % complex steerable pyramid


%% Sliding window operations
m_support = opt.pooling.m_support;
n_support = length(m_support);
n_orientations = spyr.n_orientations;
n_scales = spyr.n_scales;
if opt.do_legacy_pooling  % only for old seedoflife
  for i=1:n_scales
    n_oe = numel(spyr.oe.b{i}{1});
    for j=1:n_support
      if j ~= n_support
        continue;  % SPEED HACK
      end

      %oe_slwin{i}{j} = zeros(m_support(j)^2, n_orientations, n_oe);
      %oe_slwin = zeros(m_support(j)^2, n_orientations, n_oe);
      k_notincircle = find(nsd.util.imdisc(m_support(j),m_support(j),[],[],max(((m_support(j)-1)/2),1),1,0) == 0);  % TESTING: circle mask 
      k_incircle = find(nsd.util.imdisc(m_support(j),m_support(j),[],[],max(((m_support(j)-1)/2),1),1,0));  % TESTING: circle mask 
      for k=1:n_orientations
        %oe_slwin = zeros(m_support(j)^2, n_oe);
        b = spyr.oe.b{i}{k};

        switch lower(opt.pooling.mode)
          case 'max'
            oe_slwin = nsd.util.im2col(b, m_support(j)); % TESTING
            oe_slwin(k_notincircle,:) = 0;          
            b_zeropad = padarray(b,[m_support(j)+1 m_support(j)+1],0,'both');
            oe_slwin_zp = nsd.util.im2col(b_zeropad, m_support(j)); % TESTING
            oe_slwin_zp(k_notincircle,:) = 0;
            max_oe_slwin{i}{j}(:,k,:) = max(oe_slwin,[],1);  % max band over circular support
            pooled_oe_slwin_zp{i}{j}(:,k,:) = max(oe_slwin_zp,[],1);  % max band over circular support
            k_zeropad = 0;

          case 'sum'
            %h = [ones(m_support(j),1)];           
            %h = [1 2 3 4 3 2 1]; h = h ./ sum(h);          
            %h = [0 1 2 3 2 1 0]; h = h ./ sum(h);          
            h = nsd.util.gaussian((-3:1:3)./2.5); h = h ./ sum(h);
            max_oe_slwin{i}{j}(1,k,:) = nsd.util.vectorize(conv2(h,h,b,'same'));
            pooled_oe_slwin_zp{i}{j}(1,k,:) = nsd.util.vectorize(conv2(h,h,padarray(b,[1 1],0),'full'));
            k_zeropad = floor(length(h)/2) + 1;

            %           b_pool = conv2(h,h,b,'valid');
            %           max_oe_slwin{i}{j}(1,k,:) = nsd.util.vectorize(padarray(b_pool,[3 3],'replicate'));
            %           pooled_oe_slwin_zp{i}{j}(1,k,:) = nsd.util.vectorize(padarray(b_pool,[7 7],'replicate'));

            %oe_slwin = nsd.util.im2col(b, m_support(j)); % TESTING
            %oe_slwin(k_notincircle,:) = 0;          
            %b_zeropad = padarray(b,[m_support(j)+1 m_support(j)+1],0,'both');
            %oe_slwin_zp = nsd.util.im2col(b_zeropad, m_support(j)); % TESTING
            %oe_slwin_zp(k_notincircle,:) = 0;
            %max_oe_slwin{i}{j}(:,k,:) = sum(oe_slwin,1);  % SUM
            %pooled_oe_slwin_zp{i}{j}(:,k,:) = sum(oe_slwin_zp,1);  % SUM (zeropadded)
            %k_zeropad = m_support(j)+1;
          otherwise
            error('invalid pooling mode');
        end
        oe{i}{j}(:,k,:) = spyr.oe.b{i}{k}(:)';  % TESTING: zeroth scale      
      end
    end
  end
else
  % unused
  max_oe_slwin = [];
  pooled_oe_slwin_zp = [];
  k_zeropad = 0;
  oe = [];
end

m_scaledsupport = m_support;
m_supportscale = ones(size(m_support));
for k=2:n_scales
  m_scaledsupport = [m_scaledsupport 2*m_scaledsupport(end)+1];
  m_supportscale = [m_supportscale k];
end


%% Pooled orientation energy (new seedoflife)
oesize = size(spyr.oe.b{1}{1});  % even power of two
for i=1:n_scales
  for j=1:n_orientations
    
    switch lower(opt.pooling.mode)
      case 'sum'
        %oe_marginal(j,i,:,:) = imresize(spyr.oe.b{i}{j},oesize);

        h = nsd.util.gaussian((-3:1:3)./2.5); h = h ./ sum(h);
        b = padarray(conv2(h,h,spyr.oe.b{i}{j},'valid'),[3 3],'replicate','both');

        %h = nsd.util.gaussian((-7:1:7)./2.5); h = h ./ sum(h);  % TESTING (worse on vgg-affine)    
        %b = padarray(conv2(h,h,spyr.oe.b{i}{j},'valid'),[7 7],'replicate','both');

        oe_pooled(j,i,:,:) = imresize(b,oesize); % default 'bicubic' (nearest is worse on vgg-affine)
        %oe_marginal(j,i,:,:) = imresize(conv2(h,h,spyr.oe.b{i}{j},'same'),oesize); % default 'bicubic'
      
      case 'max'
        b = spyr.oe.b{i}{j};
        k_notincircle = find(nsd.util.imdisc(m_support(3),m_support(3),[],[],max(((m_support(3)-1)/2),1),1,0) == 0);  % TESTING: circle mask
        oe_slwin = nsd.util.im2col(b, m_support(3)); % TESTING
        oe_slwin(k_notincircle,:) = 0;
        oe_pooled(j,i,:,:) = imresize(reshape(squeeze(max(oe_slwin,[],1)),size(b)),oesize);  % max band over circular support

      otherwise
        error('invalid pooling mode ''%s''', opt.pooling.mode);
    end
  end
end


%% Features
f.img = img;
f.imgrey = imgrey;
f.imorig = imorig;
f.imscale = size(imgrey,1)./size(imorig,1);
f.spyr = spyr;
f.center = size(imgrey)./2;
f.max_oe_slwin = max_oe_slwin;
f.pooled_oe_slwin_zp = pooled_oe_slwin_zp;
f.oe = oe;
f.oe_pooled = oe_pooled;
%f.oe_slwin = oe_slwin;
f.m_support = opt.pooling.m_support;
f.m_scales = 1:n_scales;
f.m_scaledsupport = m_scaledsupport;
f.n_scaledsupport = length(m_scaledsupport);
f.m_supportscale = m_supportscale;
f.n_scales = n_scales;
f.n_support = length(opt.pooling.m_support);
f.n_orientations = n_orientations;
f.k_zeropad = k_zeropad;
f.opt = opt;
