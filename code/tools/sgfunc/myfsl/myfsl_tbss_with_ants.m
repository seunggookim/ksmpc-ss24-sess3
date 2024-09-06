function JOB = myfsl_tbss_with_ants(JOB)
% JOB = myfsl_tbss_with_ants (JOB)
% 
%
% (cc) 2015, sgKIM.  solleo@gmail.com  https://ggooo.wordpress.com

if ~isfield(JOB,'moving_suffix')
JOB.moving_suffix='';
end
if ~isfield(JOB,'fixed')
JOB.fixed = 'fa58';
end
if ~isfield(JOB,'dir_out')
JOB.dir_out = [pwd,'/ants+tbss/'];
[~,~] = mkdir(JOB.dir_out);
end

%% 1. ANTs for normalization
Deformed = myants_warp(JOB.moving_prefix, JOB.subjID, JOB.moving_suffix, JOB.fixed, JOB.dir_out);

%% 2. tbss_3_postreg: skeletonizing
[~,~]=mkdir(JOB.dir_out,'/FA');
% move all files to "FA"

% and what's inside ".msf"? final cost value?

!echo "finding best target"
!n=`$FSLDIR/bin/imglob *_FA.nii* *_FA.img* *_FA.hdr* | wc -w`
!for f in `$FSLDIR/bin/imglob *_FA.nii* *_FA.img* *_FA.hdr*` ; do
! meanscore=0
! medianscore=0
! for g in `$FSLDIR/bin/imglob *_FA.nii* *_FA.img* *_FA.hdr*` ; do
!   thismeanscore=`cat ${g}_to_${f}_warp.msf | awk '{print $1}'`
!   thismedianscore=`cat ${g}_to_${f}_warp.msf | awk '{print $2}'`
!   meanscore=`echo "10 k $meanscore $thismeanscore + p" | dc -`
!   medianscore=`echo "10 k $medianscore $thismedianscore + p" | dc -`
! done
! meanscore=`echo "10 k $meanscore $n / p" | dc -`
! medianscore=`echo "10 k $medianscore $n / p" | dc -`
! echo "$f $meanscore $medianscore"
! echo "$f $meanscore $medianscore" >> all.msf
!done
!best=`cat all.msf | sort -k 2 -n | head -n 1 | awk '{print $1}'`
!echo "best target is $best - now registering this to standard space"
!$FSLDIR/bin/imcp $best target
!echo $best > best.msf
!mkdir -p ../stats
unix('tbss_3_postreg -S');

%% 
[~,~]=mkdir([JOB.dir_out,'/stats']);
cd ([JOB.dir_out,'/stats']);
unix('tbss_4_prestats 0.2');

% %% GLM
% for k=1%:4
%   meas = LOCALMEAS{k};
%   % get links
%   dir_glm = [dir_output,'/glm/',meas,'_',JOB.glm_prefix,fsss_model_desc(JOB.model)];
%   [~,~] = mkdir(dir_glm);
%   cd(dir_glm);
%   fnames=dir([dir_output,'/stats/*',meas,'*']);
%   for j=1:numel(fnames)
%     if ~fnames(j).isdir
%       unix(['ln -sf ',dir_output,'/stats/',fnames(j).name,' ',dir_glm,'/.']);
%     end
%   end
%   unix(['ln -sf ',dir_output,'/stats/mean_FA_skeleton_mask.nii.gz ',dir_glm,'/.']);
%   
%   % get design and contrasts
%   NR=size(JOB.model,2); % # of regressors
%   varnames = char(JOB.model);
%   if ~isfield(JOB,'contrasts')
%     if isfield(JOB,'cidx')
%       cidx = JOB.cidx;
%     end
%   else
%     cidx = find(JOB.contrasts(1,:));
%   end
%   JOB.contrasts = [zeros(1,cidx-1),  1, zeros(1,NR-cidx); zeros(1,cidx-1), -1, zeros(1,NR-cidx)];
%   NumLevels = numel(unique(JOB.model(cidx)));
%   if NumLevels==2 && ~isfield(JOB,'grouplabel')
%     JOB.grouplabel={['Non-',varnames{cidx}],varnames{cidx}};
%   end
%   M = double(JOB.model);
%   fname_base = varnames{cidx};
%   myfsl_designmatrix(fname_base, M, NumLevels, JOB.grouplabel, JOB.contrasts);
%   
%   % now run randomization
%   if ~isfield(JOB,'NumPerm')
%     NumPerm=10000;
%   else
%     NumPerm=JOB.NumPerm;
%   end
%   unix([fsl_env,'cd ',dir_glm,'; randomise_parallel -i all_',meas,'_skeletonised ',...
%     ' -o ',fname_base,' -m mean_FA_skeleton_mask ', ...
%     ' -d ',fname_base,'.mat -t ',fname_base,'.con -n ',num2str(NumPerm),' --T2 -V --uncorrp'])
%   %' -d ',fname_base,'.mat -t ',fname_base,'.con -n ',num2str(NumPerm),' --T2 -V'])
% end

JOB = myfsl_randomise(JOB)



end
