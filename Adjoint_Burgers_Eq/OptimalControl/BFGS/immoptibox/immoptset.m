function  opts = immoptset(mlfun, varargin)
%IMMOPTSET  Set options for the immoptibox function MLFUN
% Call    opts = immoptset(mlfun, p1,v1, p2,v2, ...)
%
% Input parameters
% mlfun     :  String with the name of the function or a handle to it.
% p1,p2,... :  Strings with option names.
% v1,v2,... :  Real numbers with the values of the options.
%
% Output parameter
% opts :  Struct with the option names and their values.  Options 
%         that do not appear as input are assigned their default value.

% Version 10.10.19.  cv(a)imm.dtu.dk and hbn(a)imm.dtu.dk


% Define default options of respective functions:
dampnewton = struct('tau',1e-3, 'tolg',1e-4, 'tolx',1e-8, 'maxeval',100);
linesearch = struct('choice',1, 'cp1',1e-3, 'cp2',0.99, 'maxeval',10, ...
  'amax',10);
ucminf = struct('Delta',1, 'tolg',1e-4, 'tolx',1e-8, 'maxeval',100);
dogleg = struct('Delta',0, 'tolg',1e-4, 'tolx',1e-8, 'tolr',1e-6, ...
  'maxeval',100);                % NB Default Delta depends on x
marquardt = struct('tau',1e-3, 'tolg',1e-4, 'tolx',1e-8, 'maxeval',100);
smarquardt = struct('tau',1e-3, 'tolg',1e-4, 'tolx',1e-8, ...
  'maxeval',0, 'relstep',1e-7);  % NB Default maxeval depends on n
nonlinhuber = struct('gamma',0, 'Htype',0, 'tau',1e-3, 'tolg',1e-4, ...
  'tolx',1e-8, 'maxeval',100);
mexpfit = struct('const',0, 'tau',1e-3, 'tolg',1e-6, 'tolx',1e-10, ...
  'maxeval',100);
nonlinsys = struct('Delta',0, 'tolg',1e-6, 'tolx',1e-8, ...
  'maxeval',100, 'relstep',1e-6);

% List of functions:
funlist = struct('dampnewton',dampnewton, 'linesearch',linesearch, ...
  'ucminf',ucminf, 'dogleg',dogleg, 'marquardt',marquardt, ...
  'smarquardt',smarquardt, 'nonlinhuber',nonlinhuber, ...
  'mexpfit',mexpfit, 'nonlinsys',nonlinsys);
funnames = fieldnames(funlist);

% Check mlfun
if  isa(mlfun,'function_handle'),  mlfun = lower(func2str(mlfun));
elseif  ischar(mlfun),  mlfun = lower(mlfun);
else
  error('mlfun  must be a string or a function handle')
end
foundname = 0;  kfun = 0;
while  ~foundname && (kfun < length(funnames))
  kfun = kfun + 1;
  foundname = strcmp(funnames(kfun), mlfun);
end
if ~foundname
  error(['IMMOPTSET is not prepared for ', mlfun])
end

% Check existence of options and consistency of varargin                
opts = funlist.(mlfun);  optnames = fieldnames(opts);
lst = length(varargin);
for i = 1 : 2 : lst-1
  if ~ischar(varargin{i})
    error(['The value ' num2str(varargin{i}) ...
      ' is not assigned to an option'])
  end
  foundname = 0; j = 0;
  while  ~foundname && j < length(optnames)
    j = j + 1;  foundname =  strcmpi(optnames{j},varargin{i});
  end
  if  ~foundname
    error(['Option ' num2str(varargin{i}) ...
      ' does not exist in function ',mlfun])
  end
  if ischar(varargin{i+1})
    error(['No value specified for option ',num2str(varargin{i})])
  end
  % Update option
  opts.(optnames{j}) = varargin{i+1};
end
if mod(lst,2)  % odd number of elements in varargin
  if  ischar(varargin{lst})
    error('No value assigned to the last option')
  else    
    error(['The last value is not assigned to an option'])
  end
end