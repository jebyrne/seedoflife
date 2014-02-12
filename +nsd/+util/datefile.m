function [filename] = datefile(filedir,filebase,fileext)
%--------------------------------------------------------------------------
%
% Copyright (c) 2011 Jeffrey Byrne <jebyrne@cis.upenn.edu>
% $Id: datefile.m 79 2012-07-27 14:30:30Z jebyrne $
%
%--------------------------------------------------------------------------

% Inputs
if nargin < 3
  fileext = '';
end
if nargin == 0
  filedir = tempdir;
  filebase = '';
  fileext = '';
end

% My date string to force hour zero padding
mydatestr = upper(datestr(now,'ddmmmyy'));
hh = str2num(datestr(now,'HH'));
mm = str2num(datestr(now,'MM'));
if hh < 12
  mydatestr = strcat(mydatestr,sprintf('_%02d%02dAM', hh, mm));
elseif hh == 12
  mydatestr = strcat(mydatestr,sprintf('_%02d%02dPM', hh, mm));
else
  mydatestr = strcat(mydatestr,sprintf('_%02d%02dPM', hh-12, mm));
end

% Final filename
filename = fullfile(filedir, strcat(filebase, mydatestr, fileext));

