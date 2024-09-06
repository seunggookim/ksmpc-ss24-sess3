function myfsl_lcc(fn_in, thrs,fn_out)
% myfsl_lcc(fn_in, thrs,fn_out)
%
% to get a Largest Connected Component after height-thresholding
%
% (cc) 2015, sgKIM  mailto://solleo@gmail.com  https://ggooo.wordpress.com

if ~exist('thrs','var'), thrs=0.1; end

nii = load_uns_nii(fn_in);
CC = bwconncomp(nii.img>thrs);
nc = zeros(CC.NumObjects,1);
for j=1:CC.NumObjects
nc(j)=numel(CC.PixelIdxList{j});
end
[~,idx] = max(nc); % find the largest connected component
idx_lcc = CC.PixelIdxList{idx};
img = nii.img*0;
img(idx_lcc) = nii.img(idx_lcc);
nii.img = img;
if ~exist('fn_out','var')
[a,b,c] = fileparts_gz(fn_in);
fn_out = [a,'/',b,'_LCC',num2str(thrs),c];
end
save_untouch_nii(nii, fn_out);
