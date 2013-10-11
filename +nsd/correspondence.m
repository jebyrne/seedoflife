function [ij_ref, ij_refinobs, ij_obs, w_asgn, info] = correspondence(imref, imobs, opt)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne 
%
%--------------------------------------------------------------------------

%% Options
if ~exist('opt','var') || isempty(opt)
  opt = nsd.opts().correspondence;
end


%% Nested shape descriptors
[d_ref,di_ref,fr_ref,f_ref] = nsd.descriptor(imref, opt.descriptor);
[d_obs,di_obs,fr_obs,f_obs] = nsd.descriptor(imobs, opt.descriptor);
ij_ref = fr_ref(1:2,:)';
ij_obs = fr_obs(1:2,:)';


%% Pairwise Descriptor Distance
D = nsd.distance(d_ref,di_ref,fr_ref,d_obs,di_obs,fr_obs,opt.distance.mode,opt.distance);


%% Assignment
%[ij_ref, ij_refinobs, k_asgn_ref2obs, w_asgn, info] = nsd.assignment(D, fr_ref, fr_obs, opt.assignment);
switch opt.assignment.mode    
  case 'bipartite'
    if opt.verbose, fprintf('[nsd.%s]: greedy bipartite assignment \n', mfilename); end
    [u,v] = find(D);
    [dists, perm] = sort(nonzeros(D(:)),'ascend');
    %[aIdx bIdx] = ind2sub([size(d_ref,2), size(d_obs,2)], perm(1:nnz(D)));    
    k_asgn = benchmarks.helpers.greedyBipartiteMatching(size(d_ref,2), size(d_obs,2), [u(perm) v(perm)]);  % vlbenchmarks
    w_asgn = [];
    k_valid = find(k_asgn);
    ij_ref = fr_ref(1:2,k_valid)';
    ij_refinobs = fr_obs(1:2,k_asgn(k_valid))';

  case 'greedy'
    if opt.verbose, fprintf('[nsd.%s]: greedy assignment \n', mfilename); end
    [d_asgn, k_asgn] = min(D,[],2);  % one-to-many assignment, ref -> obs
    k_valid = find(isfinite(d_asgn));
    ij_ref = fr_ref(1:2,k_valid)';
    ij_refinobs = fr_obs(1:2,k_asgn(k_valid))';
    k_asgn = [k_valid k_asgn(k_valid)];
    w_asgn = exp(-d_asgn(k_valid));
        
  case 'leftright'
    fprintf('[nsd.%s]: left-right assignment \n', mfilename);
    [d_asgn, k_asgn_ref2obs] = min(D,[],2);  % one-to-one assignment, ref -> obs
    [xx, k_asgn_obs2ref] = min(D,[],1);  % one-to-one assignment, obs -> ref 
    k_valid = find(k_asgn_obs2ref(k_asgn_ref2obs) == [1:size(D,1)]);
    ij_ref = fr_ref(1:2,k_valid)';
    ij_obs = fr_obs(1:2,:)';
    ij_refinobs = ij_obs(k_asgn_ref2obs(k_valid),:);
    w_asgn = exp(-d_asgn(k_valid));
    D = D(k_valid,:);
    k_asgn = k_asgn_ref2obs(k_valid);
    
  otherwise 
    error('undefined assignment mode ''%s''', opt.mode);
end  


%% Refinement
switch opt.refinement.mode    
  case {'none'}
    % passthrough
  otherwise
    error('undefined refinement mode ''%s''', opt.mode);    
end

%% Debugging
info.f_obs = f_obs;
info.f_ref = f_ref;
info.fr_obs = fr_obs;
info.fr_ref = fr_ref;
info.k_asgn = k_asgn;
info.w_asgn = w_asgn;
info.d_ref = d_ref;
info.di_ref = di_ref;
info.D = D;
info.d_obs = d_obs;
info.di_obs = di_obs;
return;


