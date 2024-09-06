function [ind,foundstr,ind_set,foundstr_set] = ...
  myfsl_findlabel(niftiinput, roiname, atlasindex, verbose)
% given NIFTI volume, finds a logical index array for a given ROI keyword
%
% [ind,foundstr,indset,foundstrset] = myfsl_findlabel(nii_in, roiname)
%
% (cc) 2015, 2022, dr.seunggoo.kim@gmail.com
%
% SEE ASLO: MYFSL_ATLASQUERY

if ~exist('verbose','var'), verbose = true; end
fslpath = fullfile(getenv('FSLDIR'),'data','atlases');
atlas_xml = {};
atlas_nifti = {};
%% I'm not happy with combined bilateral ROIs...
% if ~isfile(fullfile(fslpath,'HarvardOxford-Cortical-Lateralized.xml'))
%   if verbose
%   warning(['You can download a lateralized version from NeuroDebian: ',...
%     'https://neuro.debian.net/install_pkg.html?p=fsl-harvard-oxford-',...
%     'cortical-lateralized-atlas'])
%   warning('For now just use BILATERAL HO cortical atlas')
%   end
%   atlas_xml{1} = xml2struct(fullfile(fslpath,'HarvardOxford-Cortical.xml'));
%   atlas_nifti{1} = fullfile(fslpath,'HarvardOxford/HarvardOxford-cort-maxprob-thr25-1mm.nii.gz');
% else
%   if verbose
%     fprintf('Using LATERALIZED HO cortical atlas...\n')
%   end
%   atlas_xml{1} = xml2struct(fullfile(fslpath,'HarvardOxford-Cortical-Lateralized.xml'));
%   atlas_nifti{1} = fullfile(fslpath,'HarvardOxford/HarvardOxford-cortl-maxprob-thr25-1mm.nii.gz');
% end
% EVEN the HPC doesn't have the lateralized atlas...
atlas_xml{1} = xml2struct(fullfile(fslpath,'HarvardOxford-Cortical.xml'));
atlas_nifti{1} = fullfile(fslpath,'HarvardOxford/HarvardOxford-cort-maxprob-thr25-1mm.nii.gz');

%%
atlas_xml{2} = xml2struct(fullfile(fslpath,'HarvardOxford-Subcortical.xml'));
atlas_nifti{2} = fullfile(fslpath,'HarvardOxford/HarvardOxford-sub-maxprob-thr25-1mm.nii.gz');
atlas_xml{3} = xml2struct(fullfile(fslpath,'Cerebellum_MNIfnirt.xml'));
atlas_nifti{3} = fullfile(fslpath,'Cerebellum/Cerebellum-MNIfnirt-maxprob-thr25-1mm.nii.gz');
atlas_xml{4} = xml2struct(fullfile(fslpath,'JHU-tracts.xml'));
atlas_nifti{4} = fullfile(fslpath,'JHU/JHU-ICBM-tracts-maxprob-thr25-1mm.nii.gz');
atlas_xml{5} = xml2struct(fullfile(fslpath,'JHU-labels.xml'));
atlas_nifti{5} = fullfile(fslpath,'JHU/JHU-ICBM-labels-1mm.nii.gz');

% unpack XML structures:
atlas_label = [];
for iatl = 1:5
  for iroi = 1:numel(atlas_xml{iatl}.atlas.data.label)
    atlas_label(iatl).index(iroi) = str2double(atlas_xml{iatl}.atlas.data.label{iroi}.Attributes.index) + 1;
    atlas_label(iatl).name{iroi} = atlas_xml{iatl}.atlas.data.label{iroi}.Text;
  end
end

% % decompress GZIP files:
% TempDir = tempname;
% c = onCleanup(@() rmdir(TempDir, 's'));
% for iatl = 1:5
%   gunzip(atlas_nifti{iatl}, TempDir)
%   [~,Fn,~] = fileparts(atlas_nifti{iatl});
%   atlas_nifti{iatl} = fullfile(TempDir, Fn);
% end

%% CHECK input
if ~exist('roiname','var')
  warning('No ROINAME given. ROINAME can be:')
  disp([atlas_label.name]')
  return
end

%% FIND atlases
if ~exist('atlasindex','var') || isempty(atlasindex)
  is_in_this_atlas = false(1,5);
  for iatl = 1:5
    is_in_this_atlas(iatl) = any(contains(atlas_label(iatl).name, roiname));
  end
  if sum(is_in_this_atlas) > 1
    if verbose
      warning('ROINAME=%s: Multiple atlases found!', roiname)
    end
  elseif sum(is_in_this_atlas) == 0
    error('ROINAME=%s: Not found in any atlases', roiname)
  end
  atlasindex = find(is_in_this_atlas);
end

%% FIND the label
info = niftiinfogz(niftiinput);
[i,j,k] = ndgrid(1:info.ImageSize(1),1:info.ImageSize(2),1:info.ImageSize(3));
xyz = ([i(:) j(:) k(:) 2*ones(numel(k(:)),1)]-1)*info.Transform.T;
ind_set = []; foundstr_set = [];
for iatl = atlasindex
  if verbose, tic; end
  index = find(contains(atlas_label(iatl).name, roiname));
  foundstr = atlas_label(iatl).name(index);
  ind = false(prod(info.ImageSize(1:3)),1);
  atlas_info = niftiinfogz(atlas_nifti{iatl});
  atlas_vol = niftireadgz(atlas_nifti{iatl});
  ijk_atlas = round(xyz*inv(atlas_info.Transform.T) + 1); % effectively NN
  for idx = index
    ijk_this = find3(atlas_vol == idx);
    ind = ind | ismember(ijk_atlas(:,1:3), ijk_this,'rows');
    % ind = ind | arrayfun(@(x) atlas_vol(ijk_atlas(x,1),ijk_atlas(x,2), ijk_atlas(x,3)) == idx, 1:size(ijk_atlas,1));
  end
  ind = reshape(ind, info.ImageSize(1:3));
  if verbose, toc, end % 3.74 sec... for 1-mm whole-brain (7M voxels)
  ind_set = cat(4,ind_set,ind);
  foundstr_set = [foundstr_set foundstr];
end

% SORT the set by # of voxels
[~,idx] = sort(squeeze(sum(sum(sum(ind_set,1),2),3)),'descend');
ind_set = ind_set(:,:,:,idx);
foundstr_set = foundstr_set(idx);

% %% delete tempfiles:
% for iatl = 1:5
%   if ~contains(atlas_nifti{iatl},getenv('FSLDIR'))
%     delete(atlas_nifti{iatl})
%   end
% end


