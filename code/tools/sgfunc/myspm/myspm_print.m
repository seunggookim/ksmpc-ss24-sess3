function myspm_print (fig_fname, fig_handle, fig_type)

error('IT DOES NOT WORK AS YOU THINK!!!!');

if ~exist('fig_handle','var')
fig_handel=gcf;
end
if ~exist('fig_type','var')
fig_type='jpg';
end

matlabbatch=cell(1);
matlabbatch{1}.spm.util.print.fname = fig_fname;
matlabbatch{1}.spm.util.print.fig.figname = fig_handle;
matlabbatch{1}.spm.util.print.opts = fig_type;
spm_jobman('run',matlabbatch)
end
