function fname = getfullpath(fname)

if ispc
  if ~strcmp(fname(2),':')
    fname = fullfile(pwd, fname);
  end
else
  if ~strcmp(fname(1),'/')
    fname = fullfile(pwd, fname);
  elseif strcmp(fname(1),'~')
    fname = fullfile(getenv('HOME'), fname(2:end));
  end
end

end