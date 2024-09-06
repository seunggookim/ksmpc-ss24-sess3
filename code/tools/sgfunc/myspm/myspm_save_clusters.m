xSPM = evalin('base','xSPM;');
XYZ  = xSPM.XYZ;

switch lower(action)
case 'thresh'
Z = xSPM.Z;

case 'binary'
Z = ones(size(xSPM.Z));

case 'n-ary'
Z       = spm_clusters(XYZ);
num     = max(Z);
[n, ni] = sort(histc(Z,1:num), 2, 'descend');
n       = size(ni);
n(ni)   = 1:num;
Z       = n(Z);

case 'current'
[xyzmm,i] = spm_XYZreg('NearestXYZ',...
spm_results_ui('GetCoords'),xSPM.XYZmm);
spm_results_ui('SetCoords',xSPM.XYZmm(:,i));

A   = spm_clusters(XYZ);
j   = find(A == A(i));
Z   = ones(1,numel(j));
XYZ = xSPM.XYZ(:,j);

otherwise
error('Unknown action.');
end

spm_write_filtered(Z, XYZ, xSPM.DIM, xSPM.M,...
sprintf('SPM{%c}-filtered: u = %5.3f, k = %d',xSPM.STAT,xSPM.u,xSPM.k));
