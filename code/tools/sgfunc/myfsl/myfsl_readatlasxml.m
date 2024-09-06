function atlas = myfsl_readatlasxml (fn_xml)
% atlas = myfsl_readatlasxml (fn_xml)
% (cc) 2019, sgKIM.

xml = xml2struct(fn_xml);
data = xml.Children(ismember({xml.Children.Name},'data')).Children;
labels = data(ismember({data.Name},'label')); % ignoring "#text"
atlas = [];
for i = 1:numel(labels)
  [~,idx] = ismember({labels(i).Attributes.Name},{'index','x','y','z'});  
  atlas(i).index = str2double(labels(i).Attributes(idx(1)).Value);
  atlas(i).index = atlas(i).index + 1; % for love of god it's zero-based..
  % no idea what these coordinates are. definitely not MNI(mm)
  % Frontal pole (48, 94, 94)? It's anterior bound is [10,74,5] mm
%   atlas(i).ijk0 = [str2double(labels(i).Attributes(idx(2)).Value), ...
%     str2double(labels(i).Attributes(idx(3)).Value), ...
%     str2double(labels(i).Attributes(idx(3)).Value)];
  atlas(i).label = labels(i).Children.Data;
end
