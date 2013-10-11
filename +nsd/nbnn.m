function [yhat,what,H] = nbnn(X, y, Z, n_nbrs, X_kdtree)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne <jebyrne@cis.upenn.edu>
%
%--------------------------------------------------------------------------

%% Exact nearest neighbor 
if ~exist(X_kdtree,'var')
  fprintf('[nsd.%s]: building kdtree\n', mfilename);
  X_kdtree = vl_kdtreebuild(X);  % deps/vlfeat
end


%% Naive bayes nearest neighbor
n_classes = length(unique(y));
n_desc = size(Z,2);
[knn, d] = vl_kdtreequery(X_kdtree, X, Z, 'NumNeighbors', n_nbrs); % deps/vlfeat
H = zeros(n_classes,n_nbrs);
for k=1:n_desc
  A = sparse(y(knn(:,k),1),[1:n_nbrs]',exp(-d(:,k)), n_classes, n_nbrs);
  h = max(A,[],2);  % max over classes
  hn = h ./ sum(h);  % class likelihood normalization \sum_k p(Xi|Y=Yk) = 1
  H(:,k) = h.*hn;  
end
[what,yhat] = sort(sum(H,2),'descend');


%% References
% Boiman, Schechman, Irani, "In defense of nearest neighbor based image classification", CVPR 2008



%% Debugging
% figure(1); imagesc(img); colormap(gray); axis image;
% hold on; plot(ij_desc(k,2),ij_desc(k,1),'r.','MarkerSize',50); hold off;
% figure(2); imagesc(nsd.preprocess(trainlist{y(knn(1,k),1)}{y(knn(1,k),2)},opt.imgsize,0)); colormap(gray); axis image;
% ij_anno = flipud(trainanno{y(knn(1,k),1)}{y(knn(1,k),2)});
% ij = [(64/6)*X(end-1,knn(1,k)),(64/6)*X(end,knn(1,k))]' + ij_anno;
% hold on; plot(ij(2),ij(1),'r.','MarkerSize',50); hold off;
% 
% 
% % best descriptor match in class
% if sum(A(i,:)) == 0,
%   continue;
% end
% [xx,k_bestnbrinclass] = max(A(i,:));
% ij_anno = flipud(trainanno{y(knn(k_bestnbrinclass,k),1)}{y(knn(k_bestnbrinclass,k),2)});
% ij = [(64/6)*X(end-1,knn(k_bestnbrinclass,k)),(64/6)*X(end,knn(k_bestnbrinclass,k))]' + ij_anno;
% figure(3); imagesc(nsd.preprocess(trainlist{y(knn(k_bestnbrinclass,k),1)}{y(knn(k_bestnbrinclass,k),2)},opt.imgsize,0)); colormap(gray); axis image;
% hold on; plot(ij(2),ij(1),'r.','MarkerSize',50); hold off; drawnow;
% 
% pause(0.02)
