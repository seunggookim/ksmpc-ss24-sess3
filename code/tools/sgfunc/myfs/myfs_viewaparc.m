function H = myfs_viewaparc(surfs, aparcdata, cfg)
% H = myfs_viewaparc(surfs, aparcdata, cfg)

if ~exist('cfg','var'), cfg = []; end
if ~isfield(cfg,'basesurf')
  fields = fieldnames(surfs);
  idx = find(ismember(fields, {'semiinflated','inflated','white','pial'}));
  if isempty(idx)
    error('No BASESURF recognized/defined')
  end
  cfg.basesurf = fields{idx(1)};
end

maps = cellfun(@(x) x.vertices(:,1)*0, surfs.(cfg.basesurf), 'uni',0);
for ihemi = 1:2
  vertlabels = surfs.aparc{ihemi}.label;
  for iroi = 1:numel(aparcdata{ihemi}.names)
    strnames = upper(cellfun(@(x) x(3:end), ...
      strrep(surfs.aparc{ihemi}.cot.struct_names,'_ROI',''), 'uni',0));    
    idx = find(ismember(strnames, upper(aparcdata{ihemi}.names{iroi})));
    if isempty(idx)
      error('structureName=%s NOT FOUND!', aparcdata{ihemi}.names{iroi})
    elseif numel(idx) > 1
      error('structureName=%s MULTIPLE INSTANCES FOUND!', ...
        aparcdata{ihemi}.names{iroi})
    end
    thislabel = surfs.aparc{ihemi}.cot.table(idx,5);
    maps{ihemi}(thislabel == vertlabels) = aparcdata{ihemi}.data(iroi);
  end
end
H = myfs_viewsurf(surfs, maps, cfg);


end