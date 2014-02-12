function [] = reprint(varargin)
%--------------------------------------------------------------------------
%
% Copyright (c) 2009-2011 Jeffrey Byrne
% $Id: reprint.m 79 2012-07-27 14:30:30Z jebyrne $
%
%--------------------------------------------------------------------------
persistent str_last;
persistent n_print;


%% Control
if isempty(varargin)
  str_last = []; return;
elseif length(varargin) == 1 && strcmp(varargin{1},'\n')
  str_last = []; fprintf('\n'); return;
end
if isempty(n_print)
  n_print = 0;
end


%% Print modulo N
%n_print = n_print + 1;
%if mod(n_print,10) ~= 0
%  return;
%else
%  n_print = 0;
%end
  

%% Print string
str = sprintf(varargin{:});
if ~isempty(str_last)
  % Common prefix?
  str_suffix = str;
  str_lastsuffix = str_last;
  for k=1:length(str_last)
    if str_last(k) == str(k)
      str_suffix = str(k+1:end);
      str_lastsuffix = str_last(k+1:end);
    else
      break;
    end
  end
  
  % Backspace and print
  fprintf(1,repmat('\b',1,length(str_lastsuffix)));  % suffix backspace
  str_escaped = strrep(str_suffix,'\','\\');  % escape backslash
  str_escaped = strrep(str_escaped,'%','%%');  % escape percent
  fprintf(str_escaped);  % print me
else
  % Print
  str_escaped = strrep(str,'\','\\');  % escape backslash
  str_escaped = strrep(str_escaped,'%','%%');  % escape percent
  fprintf(str_escaped); 
end

  
%% Store last string for next backspace if no trailing newline
if strcmp(str(end-1:end),'\n')
  % Clear last string to end reprint and final newline
  str_last = []; 
else
  str_last = str;
end

