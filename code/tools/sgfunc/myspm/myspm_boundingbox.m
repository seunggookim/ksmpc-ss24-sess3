function myspm_boundingbox(fn_img, arg2, fn_out)
% myspm_boundingbox(fn_img, thrs, [fn_out])
%
% fn_img : 3D/4D image filename (realigned & normalized)
% thrs   : a numeric threshold or 'nz','nn','fv'
%          OR used-defined bounding box [2x3] in mm
%
% SEE ALSO: SPM_GET_BBOX, SPM_GET_MATDIM
%
% NOTE1: THIS WORKS ONLY FOR IMAGES IN "CANONICAL" ORIENTATION (MNI).
% NOTE2: BOUNDING BOX IS JUST COMPUTED FOR THE FIRST VOLUME (3D OR 4D).
% THUS, ONLY USE IT FOR REALIGNED AND NORMALIZED IMAGES.
%
% (cc) 2020, sgKIM, solleo@gmail.com

%%
V = spm_vol(fn_img);
if strcmp(arg2,'canon')
  arg2 = [-78 -112 -70; 78 76 85];
end

% Get tightest bounding-box for a given threshold:
if ischar(arg2) || numel(arg2)==1
  [bb,vx] = spm_get_bbox(V, arg2);
  % Give one more voxels:
  bb = bb + [-abs(vx); abs(vx)];
else
  % Adjust 'correct' bounding box coordinates:
  P = spm_imatrix(V(1).mat);
  vx = P(7:9);
  bb0 = arg2;
  lose = @(x) ceil(abs(x)).*sign(x); % round away-from-zero (opposite of FIX)
  bb = V(1).mat*lose(spm_inv(V(1).mat)*[bb0'; 1 1]);
  bb = bb(1:3,:)';
end

% Find .mat and .dim to cover a given bounding-box and voxel size.
[mat, dim] = spm_get_matdim(V, vx, bb); % "canonical" orientation only

% All voxels' coordinates in the original volume:
[i,j,k] = ind2sub(V(1).dim, 1:prod(V(1).dim));
xyz = V(1).mat*[i; j; k; i*0+1];
xyz = xyz(1:3,:)';
ind = ...
  (bb(1,1)<=xyz(:,1)) & (xyz(:,1)<=bb(2,1)) & ...
  (bb(1,2)<=xyz(:,2)) & (xyz(:,2)<=bb(2,2)) & ...
  (bb(1,3)<=xyz(:,3)) & (xyz(:,3)<=bb(2,3));
clear xyz i j k
si = sum(ind);
pd = prod(dim);
assert(si==pd, ...
  'Inconsistent dimensions: prod(dim)=%i vs. sum(ind)=%i',pd, si)
clear si pd

%% Write image
fname_out = spm_file(V(1).fname, 'prefix','b');
if isfile(fname_out)
  delete(fname_out)
end
for i = 1:numel(V)
  Vo            = V(i);
  Vo.fname      = fname_out;
  Vo.mat        = mat;
  Vo.dim(1:3)   = dim;
  Vo.pinfo      = V(i).pinfo;
  Vo.descrip    = 'myspm - cropped to a bounding box';
  img = spm_read_vols(V(i));
  spm_write_vol(Vo, reshape(img(ind), dim));
end
if exist('fn_out','var')
  movefile(fname_out, fn_out)
end
%%
end

function TEST()
myspm_boundingbox('test3.nii',0)
myspm_boundingbox('test3.nii','nz')
myspm_boundingbox('test3.nii',[-78 -112 -70; 78 76 85])
myspm_boundingbox('test3.nii','canon')
end