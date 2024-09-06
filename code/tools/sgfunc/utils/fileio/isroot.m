function TF = isroot(Str)
% TF = isroot(Str)

Str(isspace(Str)) = [];
TF = strcmp(Str, filesep);
end