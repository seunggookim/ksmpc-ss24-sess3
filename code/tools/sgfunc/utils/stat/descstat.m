function Str = descstat(X)
Str = sprintf('min=%.4f, max=%.4f, mean=%.4f, median=%.4f, mode=%.4f, std=%.4f', ...
  min(X), max(X), mean(X), median(X), mode(X), std(X));
end