function [M] = opencv2matlab(yamlfile)
%--------------------------------------------------------------------------
%
% Copyright (c) 2013 Jeffrey Byrne
% $Id$
%
%--------------------------------------------------------------------------

% Header
n_rows = nan;
n_cols = nan;
fid = fopen(yamlfile);
for k=1:6  
  str = strtrim(fgetl(fid));
  if strfind(str,'rows:')
    n_rows = sscanf(str, 'rows: %d');
  elseif strfind(str,'cols:')
    n_cols = sscanf(str, 'cols: %d');
  elseif strfind(str,'data:')
    [xx,substr] = strtok(str,'[');    
    data = sscanf(substr(2:end), '%f,');    
    break;
  end
end
fclose(fid);

% Valid file?
if ~isnan(n_rows) && ~isnan(n_cols)
  % Remaining data
  [A] = importdata(yamlfile,',',k);
  
  % Convert to matlab datatype  
  X = [data(:); nsd.util.vectorize(A.data')];
  X(isnan(X))=[];  % final elements
  M = reshape(X, [n_cols, n_rows])';
else
  error('invalid yaml file');
end





