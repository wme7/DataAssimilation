function  opts = checkopts(fun, opts)
%CHECKOPTS  Replace illegal values by default values.

% Version 10.09.24.  cv@imm.dtu.dk  and  hbn@imm.dtu.dk

% Get default values in vector notation
dopts = immoptset(fun);   dfields = fieldnames(dopts);  
nopts = length(dfields);  dvopts = zeros(1,nopts);
for  i = 1 : nopts
  dvopts(i) = getfield(dopts,dfields{i});
end

if  isempty(opts)  % return default values
  opts = dvopts;
  
else  % check given options
  if  isstruct(opts)  % convert to array
    vopts = zeros(1,nopts);
    fields = fieldnames(opts);  
    for  i = 1 : length(fields)
      fnd = 0;  j = 0;
      while  ~fnd && j < length(dfields)
        j = j + 1;  fnd =  strcmpi(dfields{j},fields{i});
      end
      if  fnd,  vopts(j) = opts.(fields{i}); end
    end
  else,  vopts = opts; end
  
  if  strcmpi('linesearch', fun)  % special treatment
    if  vopts(1) == 0  % exact line search
      dvopts = [0 1e-3 1e-3 10 10];
    end
  end
  
  % Compare with default
  lo = length(vopts);
  for  i = 1 : min(nopts,lo)
    oi = vopts(i);
    if  isreal(oi) & ~isinf(oi) & ~isnan(oi) & oi > 0
      dvopts(i) = oi;
    end
  end
  if  lo > nopts,  dvopts = [dvopts 1]; end % for linesearch purpose
  opts = dvopts;
end