function Files = findfiles(query)
% Files = findfiles(query)
Files = dir(query);
if numel(Files)
  Files = strcat({Files.folder},'/',{Files.name})';
end

end
