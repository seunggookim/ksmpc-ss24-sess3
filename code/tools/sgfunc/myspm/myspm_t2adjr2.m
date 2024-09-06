function JOB = myspm_t2adjr2(JOB)
% JOB = myspm_t2adjr2(JOB)
% JOB.fname_tmap 
%
% R^2 = t^2 / (t^2 + df)
% t^2 = DF*R^2 / (1-R^2)
%
% adjR^2 = R^2 - p/(n-p-1)*(1-R^2) for n samples and p parameters
% (cc) 2019, sgKIM.

%{
REF: 
 http://people.duke.edu/~njs28/spmdatastructure.htm
 https://sscc.nimh.nih.gov/sscc/gangc/tr.html
 http://thestatsgeek.com/2013/10/28/r-squared-and-adjusted-r-squared/
%}


if ~isstruct(JOB)
  fn_tmap = JOB;
else
  fn_tmap = JOB.fname_tmap;
end
[p1,f1,e1] = myfileparts(fn_tmap);
load ([p1,'/SPM.mat'],'SPM');
df = SPM.xX.erdf; % effective residual degrees of freedom
t = load_untouch_nii(fn_tmap);

r2 = t;
r2.img = t.img.^2 ./ (t.img.^2 + df);
fn_r2map=[p1,'/spmR2',f1(5:end),e1];
save_untouch_nii(r2, fn_r2map);

r2adj = r2;
[~,p] = size(SPM.xX.X);
r2adj.img = r2.img.^2 - p/df*(1-r2.img.^2); % effective df because of autocorrelation
fn_r2adjmap=[p1,'/spmR2adj',f1(5:end),e1];
save_untouch_nii(r2adj, fn_r2adjmap);
end