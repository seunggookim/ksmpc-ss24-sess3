function M = myspm_cat_loadtemplatesurf(is32k, ismerged)
% M = myspm_cat_loadtemplatesurf(is32k, ismerged)
% is32k = false;
% ismerged = true;
% 
% if is32k
%   dn = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces_32k');
% else
%   dn = fullfile(spm('dir'),'toolbox','cat12','templates_surfaces');
% end
% 
% Pfwhm_all = {};
% Pfwhm_all{1} = (fullfile(dn,'lh.inflated.freesurfer.gii'));
% Pfwhm_all{2} = (fullfile(dn,'rh.inflated.freesurfer.gii'));
% 
% M0 = gifti({Pfwhm_all{1}, Pfwhm_all{2}});
% if ismerged
%   delete(Pfwhm_all{1}); delete(Pfwhm_all{2})
%   M.faces = [M0(1).faces; M0(2).faces+size(M0(1).vertices,1)];
%   M.vertices = [M0(1).vertices; M0(2).vertices];
% %   M.cdata = [M0(1).cdata; M0(2).cdata];
% else
%   M0 = M;
% end
end