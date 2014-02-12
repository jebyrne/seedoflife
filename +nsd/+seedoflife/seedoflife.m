function [d,di] = seedoflife(f, fr, opt)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne
%
%--------------------------------------------------------------------------


%% Inputs
if ~exist('opt','var')
  opt = nsd.opts().sol;
end


%% Descriptor Pooling
n_desc = size(fr,2);
n_bands = f.spyr.n_orientations;
n_scales = f.spyr.n_scales;
n_lobes = opt.n_lobes;  
n_lobescale = f.spyr.n_scales;  
th_lobe = 0:2*pi/(n_lobes):(2*pi);  % lobe orientation
r_lobe = opt.dr*ones(1,n_scales);  % lobe radius 
k_support = 3;  % 7x7 pooling 
fr_ij = fr(1:2,:)';  % frame position
fr_s = fr(3,:);  % frame scale
fr_r = fr(4,:);  % frame rotation
D = zeros(n_bands, n_lobes, n_lobescale, n_desc);  % descriptors
%Dm = ones(n_bands, n_lobes, n_lobescale, n_desc);  % descriptors
oesize = size(f.oe_pooled);
f_oe = reshape(f.oe_pooled, [oesize(1) oesize(2) oesize(3)*oesize(4)]); 
for i=1:n_lobes
  for j=1:n_lobescale
    r = opt.dr .* fr_s .* 2.^(j-1);
    th = th_lobe(i)-fr_r;
    [i_lobe,j_lobe] = pol2cart(th, r);  % (counterclockwise from vector I=1,J=0, pointing down, SOUTH, vertical)
    ij_lobe = [sign(i_lobe).*ceil(abs(i_lobe)); sign(j_lobe).*ceil(abs(j_lobe))];  % round away from zero (critical)
    ij_rlip = round(fr_ij'+ij_lobe); % rounded lobe interest points    
        
    % --- <FASTHACK>: vectorized border handling
    [k_invalid] = nsd.util.notinmat(oesize(3:4),ij_rlip(1,:),ij_rlip(2,:));  % invalid rounded scaled interest points?
    if ~isempty(k_invalid)
      % Closest intersection with image border of radial lobe
      p_sip = fr_ij(k_invalid,:)';  % scaled interest point center
      p_lobe = [fr_ij(k_invalid,:)+ij_lobe(:,k_invalid)']';  % scaled lobe center
      p_border = [];
      p_border(:,:,1) = nsd.util.line_line_intersection([1 1],[oesize(3) 1],p_sip,p_lobe);    % top
      p_border(:,:,2) = nsd.util.line_line_intersection([1 1],[1 oesize(4)],p_sip,p_lobe);    % left
      p_border(:,:,3) = nsd.util.line_line_intersection([oesize(3) 1],[oesize(3) oesize(4)],p_sip,p_lobe); % right
      p_border(:,:,4) = nsd.util.line_line_intersection([1 oesize(4)],[oesize(3) oesize(4)],p_sip,p_lobe); % bottom
      d_border = ndsum((repmat(p_sip,[1 1 4])-p_border).^2,1);
      p_rborder = round(p_border); % rounded
      [xx,j_asgn] = min(d_border,[],2);  % closest border intersection (top,left,right,bot)
      for k=1:length(k_invalid)
        ij_rlip(:,k_invalid(k)) = min(max(p_rborder(:,k,j_asgn(k)),1),oesize(3:4)');  % truncated rounding of lobe-border intersection
      end
    end
    [k_rlip] = nsd.util.sub2ind(oesize(3:4), ij_rlip(1,:), ij_rlip(2,:));  % rounded lobe interest points (all guaranteed valid from border handling)    
    d = f_oe(:,:,k_rlip);         
    % --- </FASTHACK>: vectorized border handling
    
    % Subband scale interpolation (vectorized)
    si = log2(fr_s);  % scale interpolation
    w = nsd.util.vectorize(repmat(ceil(abs(si)) - abs(si),n_bands,1));
    jmin = min(max(j-1+sign(si).*floor(abs(si)),1),n_scales);  % minus 1 is critical!
    jmax = min(max(j-1+sign(si).*ceil(abs(si)),1),n_scales);
    k_maxscale = sub2ind(size(d),nsd.util.vectorize(repmat([1:n_bands]',1,n_desc)),...
      nsd.util.vectorize(repmat(jmax,n_bands,1)),nsd.util.vectorize(repmat(1:n_desc,n_bands,1)));
    k_minscale = sub2ind(size(d),nsd.util.vectorize(repmat([1:n_bands]',1,n_desc)),...
      nsd.util.vectorize(repmat(jmin,n_bands,1)),nsd.util.vectorize(repmat(1:n_desc,n_bands,1)));
    D(:,i,j,:) = reshape((1-w).*d(k_maxscale) + w.*d(k_minscale), n_bands, n_desc);        
  end
end


%% Rotation intepolation
fr_r = fr(4,:);  % frame rotation
if any(fr_r > 0)
  V = reshape(D, [n_bands n_lobes*n_scales n_desc]);
  [X,Y] = meshgrid(1:n_lobes*n_scales, 1:n_bands);
  for k=1:n_desc
    if f.opt.spyr.do_signed_orientation
      ri = fr_r(k)*(n_bands/(2*pi));  
    else
      ri = fr_r(k)*(n_bands/(pi));  
    end
    Yi = mod((Y+ri)-1, n_bands) + 1;      
    v = [V(:,:,k); V(1,:,k)]; % with circular boundary 
    Vi = interp2(v, X, Yi);  % cannot steer because this is OE
    D(:,:,:,k) = reshape(Vi, [n_bands n_lobes n_scales]);    
  end
end


%% Cumulative weighting
for k=1:n_scales
 D(:,:,k,:) = D(:,:,k,:)*(n_scales+1-k);
end


%% Binarize - Logarithmic spiral 
if opt.do_logspiral
  Sa = (circshift(D,[0 1 1 0]));  
  Sb = (circshift(D,[0 1 -1 0]));
  %S = D - (0.5*Sa + 0.5*Sb);  % weighted combination 
  S = D - Sa;  % single log spiral
%  S(:,:,1,:) = D(:,:,1,:);  % smallest scale unnormalized (do not use circular shift back to one!)
  S(:,:,1,:) = -D(:,:,1,:)+Sb(:,:,1,:);  % smallest scale reverse normalized
%  S(:,:,end,:) = D(:,:,end,:)-Sa(:,:,end,:);  % smallest scale reverse normalized
 % Dm = (circshift(Dm,[0 1 1 0]));  % mask
 
if opt.binarize
   D = double(sign(S) > 0);
 else
   D = S;
 end
 %D = S;  warning('binarization disabled');
 %warning('logspiral disabled');
end


%% Vectorize
d = reshape(D, [n_bands*n_lobes*n_scales n_desc]);
%dm = reshape(Dm, [n_bands*n_lobes*n_scales n_desc]);


%% Descriptor information
clear di;
di.n_bands = n_bands;
di.n_desc = n_desc;
di.n_scales = n_scales;    
di.n_lobes = n_lobes;
di.r_nesting = 2.^[1:n_scales];
di.th = th;
di.dr = opt.dr;
di.k_support = k_support;
di.is_signed = f.spyr.is_signed;
di.ij = fr_ij';
di.fr = fr;
di.xy = flipud(di.ij);
di.k = nsd.util.sub2ind(size(f.imgrey),round(di.ij(1,:)'),round(di.ij(2,:)'));



