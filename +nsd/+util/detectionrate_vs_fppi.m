function [eval] = detectionrate_vs_fppi(tpfp, conf, n_obj, n_img, plotstyle)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne <jebyrne@cis.upenn.edu>
% $Id: detectionrate_vs_fppi.m 79 2012-07-27 14:30:30Z jebyrne $
%
%--------------------------------------------------------------------------

% %% Evaluation
% n_obj = length(find([y_anno.label] == 1));
% eval.n_test = length(y_anno);
% eval.label = zeros(eval.n_test,1);  % obj=1/notobj=0
% eval.tpfp = -ones(eval.n_test,1);   % true positive=1/false positive=-1
% eval.conf = zeros(eval.n_test,1);   % confidence
% for k=1:eval.n_test
%   eval.label(k) = y_anno(k).label;
%   if ((y_anno(k).label == 1) & (y_hat(k) == 1))
%     eval.tpfp(k) = 1;  % true positive
%   end
%   if y_hat(k) == 1
%     eval.conf(k) = double(1/conf(k));
%   else
%     eval.conf(k) = double(-1/conf(k));
%   end  
% end


%% Plots
if ~exist('plotstyle','var')
  plotstyle = 'b.-';
end

% Sorted confidence
[sc,si] = sort(conf,'descend'); 

% True/False positives
tp = zeros(length(tpfp),1);
tp(tpfp == 1) = 1;
tp = tp(si); % reordered
fp = zeros(length(tpfp),1);
fp(tpfp == -1) = 1;
fp = fp(si);  % reordered

% Detection rate, false positives per image
tp = cumsum(tp);
fp = cumsum(fp);
detection_rate = tp / n_obj;
fppi = fp / n_img;

% Plot me
h=figure(sum(char(mfilename)));
set(h,'NumberTitle','off','Toolbar','none','Menubar','none','Name','Evaluation');
plot(fppi,detection_rate,plotstyle);
xlabel('False positives per image');
ylabel('Detection rate');
axis([0 max(fppi) 0 1]);
grid on;

% NOTE: multiple plots can be used with hold on and legend

