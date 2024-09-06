function conn = myfs_surf2conn(surfs)
% conn = myfs_surf2conn(surfs)
%
% inspired by fieldtrip/private/triangle2connectivity.m
%
% (cc) 2022, dr.seunggoo.kim@gmail.com

if isstruct(surfs)
  surfs = surfs.(surfs.basesurf);
elseif ~iscell(surfs)
  surfs = {surfs};
end

n = 0; faces = [];
for i = 1:numel(surfs)
  faces = [faces; surfs{i}.faces + n];
  n = n + size(surfs{i}.vertices,1);
end
assert(n == max(faces(:)), '#verts =/= max(vertId)')
edg = unique([faces(:,[1 2]); faces(:,[2 3]); faces(:,[1 3])], 'rows');
k = size(edg,1);
% indexing/filling a sparse matrix takes forever:
% create a new one with indices:
conn = sparse([edg(:,1);edg(:,2)],[edg(:,2);edg(:,1)],true(2*k,1));

end