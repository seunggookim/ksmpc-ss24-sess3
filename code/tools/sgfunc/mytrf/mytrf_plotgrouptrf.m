function job = mytrf_plotgrouptrf(job)
% job = mytrf_plotgrouptrf(job)
%
% job requires:
%  .fnames
%  .suffix
%  .dn_out
%

[~,~] = mkdir(job.dn_out);
job.fn_out = [job.dn_out,'/TRFs_group_',job.suffix,'.pdf'];
if isfile(job.fn_out)
  return
end

%% Average across outer-CV:
nsub = numel(job.fnames);
mse_opt = []; mse_eva = []; r_opt = []; r_eva = []; lambda_opt = [];
w = []; b = [];
for isub = 1:nsub
  this = load(job.fnames{isub},'CV','job');
  mse_opt = cat(3,mse_opt, this.CV.mse_opt );
  r_opt = cat(3,r_opt, this.CV.r_opt);
  mse_eva = cat(3,mse_eva, this.CV.mse_eva );
  r_eva = cat(3,r_eva, this.CV.r_eva);
  lambda_opt = [lambda_opt; geomean(this.CV.lambda_opt)];
  w = cat(5, w, cat(4,this.CV.model.w));
  b = cat(4, b, cat(3,this.CV.model.b));
end
t = this.CV.model(1).t;
lambdas = this.job.lambdas;
chanlocs = this.job.chanlocs;

%% save them:
save([job.fn_out(1:end-3),'mat'], 'mse_opt', 'mse_eva', ...
  'r_opt', 'r_eva', 'lambdas', 'lambda_opt', 'chanlocs', 't', 'w', 'b');
if isfield(this.job,'W')
  rca = struct('k',this.job.rca, 'W',this.job.W, 'A',this.job.A);
  save([job.fn_out(1:end-3),'mat'], '-append','rca');
end

%% Plot them:
if exist('rca','var')
  h = mytrf_plotmdlrca(mse_opt, mse_eva, r_opt, r_eva, ...
    lambdas, lambda_opt, chanlocs, t, w, rca);
else
  if isfield(this.job,'band')
    h = mytrf_plotmdlband(mse_opt, mse_eva, r_opt, r_eva, ...
      lambdas, lambda_opt, chanlocs, t, w);
  else
    h = mytrf_plotmdl(mse_opt, mse_eva, r_opt, r_eva, ...
      lambdas, lambda_opt, chanlocs, t, w);
  end
end
printpdf(job.fn_out);
close(h)

end