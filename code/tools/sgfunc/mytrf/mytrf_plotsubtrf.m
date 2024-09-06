function job = mytrf_plotsubtrf(job)
% job = mytrf_plotsubtrf(job)
%
% job requires:
%  .dn_mdl
%  .dn_out
%

[p1,f1,~] = fileparts(job.fn_cv);
[~,~] = mkdir(p1);
job.fn_out = fullfile(p1,[f1,'.pdf']);
if isfile(job.fn_out)
  return
end

% Average across outer-CV:
nsub = numel(job.fn_results);
mse_opt = []; mse_eva = []; r_opt = []; r_eva = []; lambda_opt = [];
w = []; b = [];
for isub = 1:nsub
  this = load(job.fn_results{isub},'CV','job');
  mse_opt = cat(3,mse_opt, mean(this.CV.mse_udf) - this.CV.mse_opt );
  r_opt = cat(3,r_opt, this.CV.r_opt);
  mse_eva = cat(3,mse_eva, mean(this.CV.mse_udf) - this.CV.mse_eva );
  r_eva = cat(3,r_eva, this.CV.r_eva);
  lambda_opt = [lambda_opt, mean(this.CV.lambda_opt)];
  w = cat(4, w, this.CV.model.w);
  b = cat(4, b, this.CV.model.w);
end
t = this.CV.model(1).t;
lambdas = this.job.lambdas;
nfeats = size(w,1);

%%
mytrf_plotmdl(w,lambdas)



end