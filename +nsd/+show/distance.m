function [h] = distance(imobs,imref,ij_ref,ij_refinobs,ij_obs,d_obs,d_ref,di_obs,di_ref,fr_obs,fr_ref,f_obs,f_ref,D)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne
% $Id: distance.m 165 2013-08-02 21:39:17Z jebyrne $
%
%--------------------------------------------------------------------------
close all; drawnow;


%% Interactive matching display
while(1)  
  % Matching for interactive selection
  h1=figure(1); set(h1,'NumberTitle','off','Name','Matching');
  nsd.show.matching(imobs,imref,ij_ref,ij_refinobs,ij_obs,h1);

  % User selects bad match
  fprintf('[nsd.%s]: select interest point in reference\n', mfilename);
  figure(h1); k_ref = nsd.show.nearest_interestpoint(imref,ij_ref',subplot(1,2,2));
  fprintf('[nsd.%s]: select correctly matching interest point in observation\n', mfilename);
  figure(h1); k_obs = nsd.show.nearest_interestpoint(imobs,fr_obs(1:2,:),subplot(1,2,1));

  % Selected frames
  fprintf('[nsd.%s]: matched frames\n', mfilename);
  k_ref = find(ij_ref(k_ref,1)==ij_ref(:,1) & ij_ref(k_ref,2)==ij_ref(:,2));  % multiple orientations?
  k_obs = find(ij_obs(k_obs,1)==ij_obs(:,1) & ij_obs(k_obs,2)==ij_obs(:,2));  % multiple orientations?
  fr_ref_k = fr_ref(:,k_ref);   % selected match
  fr_obs_k = fr_obs(:,k_obs);   % selected match
  k_ref = k_ref(1);  % multiple orientations?
  k_obs = k_obs(1);  % multiple orientations?
  
  % there is a problem visualizing matches with multiple orientations and multiple scales
  % matches may look wrong, but are matched correctly at least one orientation
  
  
  % Greedy assignment
  [xx,k_asgn] = min(D,[],2);  
  
  % Matched frames
  h8=figure(8); set(h8,'NumberTitle','off','Name','Minimum (magenta) vs. Clicked (cyan) match');
  figure(h8); subplot(1,2,1); imagesc(imobs); colormap(gray); axis image; axis off;
  fr_sift = [fr_obs(2,k_asgn(k_ref)); fr_obs(1,k_asgn(k_ref)); 10*fr_obs(3,k_asgn(k_ref)); fr_obs(4,k_asgn(k_ref))];
  hold on;  h_sift = vl_plotframe(fr_sift); hold off;
  set(h_sift,'color','m','linewidth',3);
  figure(h8); subplot(1,2,2); imagesc(imref); colormap(gray); axis image; axis off;
  fr_sift = [fr_ref(2,(k_ref)); fr_ref(1,(k_ref)); 10*fr_ref(3,(k_ref)); fr_ref(4,(k_ref))];
  hold on;  h_sift = vl_plotframe(fr_sift); hold off;
  set(h_sift,'color','m','linewidth',3);
  figure(h8); subplot(1,2,1); 
  fr_sift = [fr_obs(2,k_obs); fr_obs(1,k_obs); 10*fr_obs(3,k_obs); fr_obs(4,k_obs)];
  hold on;  h_sift = vl_plotframe(fr_sift); hold off;
  set(h_sift,'color','c','linewidth',3);
  
  % Descriptor distance 
  fprintf('[nsd.%s]: distance\n', mfilename);
  h3=figure(3);  set(h3,'NumberTitle','off','Name','Reference Descriptor Pairwise Distance');
  plot(D(k_ref,:));  title('distance');
  fprintf('[nsd.%s]: minimum match distance = %f, selected match distance = %f\n', mfilename,min(D(k_ref,:)),D(k_ref,k_obs));
  
  % (Selected match - true match) descriptor relative error
  z_true = (d_ref(:,k_ref)-d_obs(:,k_obs)).^2;  % manually selected
  z_est = (d_ref(:,k_ref)-d_obs(:,k_asgn(k_ref))).^2; % greedy assignment
  Z_true = nsd.seedoflife.reshape(z_true,di_ref);
  Z_est = nsd.seedoflife.reshape(z_est,di_ref);
  h12=figure(12); set(h12,'NumberTitle','off','Name','Descriptor Relative Error (selected-true)');
  figure(h12); subplot(2,2,1); plot(1:size(d_ref,1),z_true-z_est,'b-');  title('descriptor relative error'); legend({'(clicked-min)'}); grid on;
  figure(h12); subplot(2,2,2); imagesc(ndsum(Z_true,1)-ndsum(Z_est,1));  colorbar; title('lobes vs scale');
  figure(h12); subplot(2,2,3);imagesc(ndsum(Z_true,2)-ndsum(Z_est,2));  colorbar; title('bands vs scale');
  figure(h12); subplot(2,2,4);imagesc(ndsum(Z_true,3)-ndsum(Z_est,3));  colorbar; title('bands vs lobes');
    
  % Reference -> True observation descriptor error 
  h11=figure(11); set(h11,'NumberTitle','off','Name','Descriptor Error (true)');
  figure(h11); subplot(2,2,1); plot(1:size(d_ref,1),z_true,'b-');  title('descriptor error'); grid on;
  figure(h11); subplot(2,2,2); imagesc(ndsum(Z_true,1)); colorbar; title('lobes vs scale');
  figure(h11); subplot(2,2,3); imagesc(ndsum(Z_true,2));  colorbar; title('bands vs scale');
  figure(h11); subplot(2,2,4); imagesc(ndsum(Z_true,3));  colorbar; title('bands vs lobes');

  % Reference descriptor visualization
  h10=figure(10); set(h10,'NumberTitle','off','Name','Descriptor');
  D_ref = nsd.seedoflife.reshape(d_ref(:,k_ref),di_ref);
  figure(h10); subplot(2,2,1); plot(1:size(d_ref,1),d_ref(:,k_ref),'b-');  title('descriptor'); grid on;
  figure(h10); subplot(2,2,2); imagesc(ndsum(D_ref,1)); colorbar; title('lobes vs scale');
  figure(h10); subplot(2,2,3); imagesc(ndsum(D_ref,2));  colorbar; title('bands vs scale');
  figure(h10); subplot(2,2,4); imagesc(ndsum(D_ref,3));  colorbar; title('bands vs lobes');
  
end

