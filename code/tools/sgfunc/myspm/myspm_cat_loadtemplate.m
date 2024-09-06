function M = myspm_cat_loadtemplatesurf(is32k, ismerged)
% myspm_cat_loadtemplate(is32k, ismerged)




        M0 = gifti({Pfwhm_all{1}, Pfwhm_all{2}});
        delete(Pfwhm_all{1}); delete(Pfwhm_all{2})
        M.faces = [M0(1).faces; M0(2).faces+size(M0(1).vertices,1)];
        M.vertices = [M0(1).vertices; M0(2).vertices];
        M.cdata = [M0(1).cdata; M0(2).cdata];
        
end