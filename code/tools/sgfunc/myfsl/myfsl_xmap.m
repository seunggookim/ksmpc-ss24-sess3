function myfsl_xmap (fname_vmap,flag)
%  myfsl_xmap (fname_vmap, Log_flag)
%
% normalizes and standardizes visitaion maps into cross(X)-maps
%
% (cc) 2015, sgKIM.
if nargin<1, help myfsl_xmap; end
if ~exist('flag','var'), flag='n'; end

[path1,name1,ext1] = fileparts(fname_vmap);
waytotal = load([path1,'/waytotal']);
nii = load_uns_nii(fname_vmap);
switch lower(flag)
case 'n' % n-map
nii.img = (nii.img)./(waytotal);  % just normalized by total # of tracts... but it's greater than one?!!!
save_untouch_nii(nii, [path1,'/n',name1(2:end),ext1]);
case 'c' % c-map
nii.img = (nii.img)./(waytotal);  % just normalized by total # of tracts... but it's greater than one?!!!
% here remove small values: 0.05 looks okay..
CC = bwconncomp(nii.img>0.05);
nc = zeros(CC.NumObjects,1);
for j=1:CC.NumObjects
nc(j)=numel(CC.PixelIdxList{j});
end
[~,idx] = max(nc); % find the largest connected component
idx_lcc = CC.PixelIdxList{idx};
img = nii.img*0;
img(idx_lcc) = nii.img(idx_lcc);
nii.img = img;
save_untouch_nii(nii, [path1,'/c',name1(2:end),ext1]);
case 'l' % l-map
% logarithm from power dist to somewhat more similar to a central dist (?)
nii.img = log10(nii.img)./log10(waytotal);  % normalized by total # of tracts
nii.img(isinf(nii.img))=0;
nii.img(isnan(nii.img))=0;
save_untouch_nii(nii, [path1,'/l',name1(2:end),ext1]);
end

end
