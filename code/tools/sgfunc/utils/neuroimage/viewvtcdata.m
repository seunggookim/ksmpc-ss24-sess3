function [H,cfg] = viewvtcdata(data, InfoVTC, cfg)
% viewvtcdata( data, InfoVTC, cfg )
if ~exist('cfg','var'), cfg=[]; end
map = zeros(InfoVTC.DimVTC);
map(InfoVTC.voxVTC) = data;

switch numel(InfoVTC.DimVTC)
  case 3 % volume
    % on which template?
    if ~isempty(getenv('FSL_DIR'))
      t1w = fullfile(getenv('FSL_DIR'),'data','standard','MNI152_T1_1mm.nii.gz');
    elseif ~isempty(which('spm'))
      t1w = fullfile(spm('dir'),'canonical','avg152T1.nii');
    else
      error('CANNOT FIND NEITHER FSL NOR SPM!')
    end
    map = struct('vol',map, 'vox2ras', vox2ras_1to0(InfoVTC.vox2ras1));
% map = struct('vol',map, 'vox2ras', (InfoVTC.vox2ras1));
    slices(t1w, map, cfg);
%     myfs_viewvol(t1w,map,cfg)
    
  case 2 % surface
    error('% okay on which surface???')
end

end
