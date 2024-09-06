function Str = formatpval(P, IsShort)
if not(exist('IsShort','var'))
  IsShort = 0;
end
if P < 0.001
  if IsShort == 1
    Str = sprintf('%.0i', P); % for SHORT TEXT
  elseif IsShort == 2
    Str = sprintf('< 10^{%i}', ceil(log10(P))); % for PLOT
  else
    Str = sprintf('< 10^%i', ceil(log10(P))); % for TABLE
  end
else
  Str = sprintf('%.3f', P);
end

end
