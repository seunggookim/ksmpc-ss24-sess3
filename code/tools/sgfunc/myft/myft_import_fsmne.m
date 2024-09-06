function myft_import_fsmne(subjects_dir, subj, fn_rawmeg, dn_out)
% imports Head Model and Source Spaces from MNE to use in FieldTrip
% 
% myft_import_fsmne(subjects_dir,subj, fn_rawmeg, dn_out)
%
% REF: http://www.fieldtriptoolbox.org/example/import_mne/
% NOTE: This script was made by Mikkel Vinding, based on code courtesy of Bushra Riaz Syeda
% (cc) modified by sgKIM, 2019.

%% PRESET FILENAMES (for my pipeline):
if ~exist('subjects_dir','var')
  subjects_dir = getenv('SUBJECTS_DIR');
else
  setenv('SUBJECTS_DIR',subjects_dir);
end
fn_T1 = [subjects_dir,'/',subj,'/mri/T1.mgz'];
fn_trans = [subjects_dir,'/',subj,'/meg/',subj,'-trans.fif'];
fn_src = [subjects_dir,'/',subj,'/bem/',subj,'-4098-mnec-src.fif'];
if ~exist(fn_src,'file')
  % NOTE: some of MNE-python do not work with MNE-C/MATLAB functions: use MNE-C
  unix(['mne_make_source_space --subject ',subj,...
    ' --ico -6 --src ',fn_src,' --surf lh.white:rh.white'])
end
fn_bem = [subjects_dir,'/',subj,'/bem/',subj,'-4098-1layer-bem.fif'];
fn_ftsrc = [dn_out,'/src-surf.mat'];
fn_fthead = [dn_out,'/head-surf.mat'];
fn_png = [dn_out,'/meg-mri-coreg.png'];

%% READ head-to-MRI transformation from MNE
trans_orig = fiff_read_coord_trans(fn_trans);
%
% In FieldTrip every thing is in head coordinates:
trans.trans = inv(trans_orig.trans);    % now MRI-to-head transform

%% Read MRI (does not work on Win PC)
mri_orig = ft_read_mri(fn_T1,'dataformat','freesurfer_mgz');
mri_orig = ft_convert_units(mri_orig, 'cm');

% The following lines are importing MNE coregistration to FieldTrip:
trans_cm = trans_orig.trans;
trans_cm(1:3,4) = trans_cm(1:3,4)*100;  % translation: meters to cm
ttt = mri_orig.hdr.tkrvox2ras;          % This is for FS T1.mgz!
ttt(1:3,:) = ttt(1:3,:)/10;
mri_orig.transform = inv(trans_cm) * ttt;
mri_orig.coordsys = 'neuromag';         % It's now in head coordsys (RAS)
% mri_orig = ft_determine_coordsys(mri_orig, 'interactive', 'no');

% Read sensor and headpoints
sens = ft_read_sens(fn_rawmeg, 'senstype','meg');
sens = ft_convert_units(sens, 'cm');
hpnt = ft_read_headshape(fn_rawmeg);
hpnt = ft_convert_units(hpnt, 'cm');

%% Read headmodel (BEM)
[head] = ft_read_vol(fn_bem);
[head] = ft_convert_units(head,'m'); % Make sure units is in meters for transform

% Transform
headmodel_pos = head.bnd.pos;
temp_vect = headmodel_pos;
temp_vect(:,4) = 1;
headmodel_pos = temp_vect*trans.trans';
head.pos = headmodel_pos(:,1:3);
head.bnd.pos = head.pos; % this is what ft_prepare_leadfield reads
[head] = ft_convert_units(head, 'cm');

%% Reading FreeSurfer Source Spaces
[src] = ft_read_headshape(fn_src, 'format', 'mne_source');

% Transform
temp_vect = src.pos;
temp_vect(:,4) = 1;
src_pos = temp_vect*trans.trans';
src_pos = src_pos(:,1:3);
src.pos = src_pos;

temp_vect = src.orig.pos;
temp_vect(:,4) = 1;
src_pos = temp_vect*trans.trans';
src_pos = src_pos(:,1:3);
src.orig.pos = src_pos;

% Save a source model
[src] = ft_convert_units(src, 'cm');
src.inside = ones(length(src.pos),1);

%% Check coregistration
figure('color','k', 'position',[2524 277 896 483]);
views = {[90 0],[180 0]};
viewnames = {'lat','front'};
for i = 1:2
  subplot(1,2,i)
  ft_determine_coordsys(mri_orig, 'interactive', 'no');
  hold on;
  ft_plot_headshape_sg(hpnt);
  ft_plot_sens_sg(sens, 'style', '*g','edgecolor','cyan');
  ft_plot_mesh(head, 'facealpha',0.25, 'edgecolor','b', 'facecolor','b', ...
    'vertexcolor','m')
  ft_plot_mesh(src, 'edgecolor','k')
  view(views{i});
  title(['MEG-MRI coregistration: ',subj,' (',viewnames{i},')'], ...
    'FontSize', 13, 'color','w')
end
drawnow
export_fig(fn_png,'-r300')
close(gcf)

%% Save stuff
disp('Saving...');
save(fn_ftsrc, 'src', '-v7.3');
save(fn_fthead, 'head', '-v7.3');
disp('DONE')
end